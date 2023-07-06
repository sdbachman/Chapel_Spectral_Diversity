use distance_mask;
use distance;
use BlockDist;
use Time;
use IO_Module;
use IO;
use AutoMath;
use LinearAlgebra;
use AllLocalesBarriers;

/* Command line arguments. */
config const in_name : string;
config const window_size : real(32);
config const dx : real(32) = 2.0;

var bs = 1;
var be = 5;

proc convolve_and_calculate(Array: [] real(32), centerPoints : ?, locL : ?, locC : ?, locR : ?, Output: [] real, Mask_Size : int, t: stopwatch) : [] {

  var first_point = centerPoints.first[1];
  var last_point = centerPoints.last[1];

  forall i in centerPoints[..,first_point] {

    var tmpLL : real = 0;
    var tmpLC : real = 0;
    var tmpLR : real = 0;
    var tmpCC : real = 0;
    var tmpCR : real = 0;
    var tmpRR : real = 0;

    calc_distance(Array, locL, locL, tmpLL, bs, be, i, first_point);
    calc_distance(Array, locL, locC, tmpLC, bs, be, i, first_point);
    calc_distance(Array, locL, locR, tmpLR, bs, be, i, first_point);
    calc_distance(Array, locC, locC, tmpCC, bs, be, i, first_point);
    calc_distance(Array, locC, locR, tmpCR, bs, be, i, first_point);
    calc_distance(Array, locR, locR, tmpRR, bs, be, i, first_point);

    Output[i,first_point] = (tmpLL + tmpCC + tmpRR + 2*(tmpLC + tmpLR + tmpCR)) / (Mask_Size**2);
    var prev = tmpCC + tmpRR + 2*tmpCR;

    for j in (first_point+1)..last_point do {

      tmpLL = 0;
      tmpLC = 0;
      tmpLR = 0;
      tmpCR = 0;
      tmpRR = 0;

      calc_distance(Array, locL, locL, tmpLL, bs, be, i, j);
      calc_distance(Array, locL, locC, tmpLC, bs, be, i, j);
      calc_distance(Array, locL, locR, tmpLR, bs, be, i, j);
      calc_distance(Array, locC, locR, tmpCR, bs, be, i, j);
      calc_distance(Array, locR, locR, tmpRR, bs, be, i, j);

      var current = tmpRR + 2*(tmpLR + tmpCR);
      Output[i,j] = (prev + current) / (Mask_Size**2);
      prev = prev + current - tmpLL - 2*(tmpLC + tmpLR);

    }
  }

  writeln("Elapsed time on ", here.name, ": ", t.elapsed(), " seconds for domain ", centerPoints);

}


proc main(args: [] string) {

  var t : stopwatch;
  t.start();

  const radius = (sqrt(window_size) / 2) : int;
  const nx = (radius / dx) : int;
  writeln("Distance circle has a radius of ", nx, " points.");

  /////////////////////////////////////////////////////////////////////////////////
  // Open the text file that contains the number of rows and columns for the image.
  const image_infile = open(in_name + ".txt", iomode.r);
  const image_reader = image_infile.reader();
  // Read the number of rows and columns in the array in from the file.
  const rows = image_reader.read(int);
  const cols = image_reader.read(int);

  const ImageSpace = {0..<rows, 0..<cols};

  // Close the file.
  image_reader.close();
  image_infile.close();
  ////////////////////////////////////////////////////////////////////////////////

  var Array : [0..<rows, 0..<cols, 1..5] real(32);

  // Read in array
  var f = open(in_name + ".bin", iomode.r);
  var r = f.reader(kind=ionative);

  for i in ImageSpace.dim[0] {
    for j in ImageSpace.dim[1] {
      for k in 1..5 {
        var tmp : real(32);
        r.readBinary(tmp);
        Array[i,j,k] = tmp;
      }
    }
  }
  r.close();


  // Create distance mask
  var (LeftMask, CenterMask, RightMask, Mask_Size) = create_distance_mask(radius, dx, nx);

  // Create Block distribution of interior of PNG
  const offset = nx; // maybe needs to be +1 to account for truncation?
  const Inner = ImageSpace.expand(-offset);
  const myTargetLocales = reshape(Locales, {1..Locales.size, 1..1});
  const D = Inner dmapped Block(Inner, targetLocales=myTargetLocales);
  var OutputArray : [D] real;


  // Create NetCDF
  var varid : int;
  var Ha = window_size / 10000 / (4.0 / pi);
  var zero_pad = 3 - (Ha : string).size;
  var nonzero : string = "%3.1dr".format(Ha);
  var Ha_string : string = (zero_pad*"0") + nonzero;

  var out_file = in_name + "_" + Ha_string + "Ha.nc";
  CreateNetCDF(out_file, ImageSpace, varid);

  writeln("Elapsed time at start of coforall loop: ", t.elapsed(), " seconds.");

  writeln("Starting coforall loop.");

  coforall loc in Locales do on loc {

    var t2 : stopwatch;
    t2.start();

    // If I put "create_distance_mask" inside this loop I need to declare local copies of these variables,
    // otherwise it seems like Chapel will have to do a ton
    // of cross-locale calls to access these variables. This seems to double the amount of time it
    // takes to run through the coforall loop for all non-head locales!

    const loc_Mask_Size = Mask_Size;

    const locArrayDomain = Array.domain;
    const locArray : [locArrayDomain] Array.eltType = Array;

    const locLeftMaskDomain = LeftMask.domain;
    const locCenterMaskDomain = CenterMask.domain;
    const locRightMaskDomain = RightMask.domain;

    allLocalesBarrier.barrier();
    convolve_and_calculate(locArray, D.localSubdomain(), locLeftMaskDomain, locCenterMaskDomain, locRightMaskDomain, OutputArray, loc_Mask_Size, t2);
    t2.stop();

    writeln("Elapsed time to calculate Q on ", here.name, ": ", t2.elapsed(), " seconds.");
  }


  writeln("Elapsed time to finish coforall loop: ", t.elapsed(), " seconds.");

  WriteOutput(out_file, OutputArray, varid);

  writeln("Elapsed time to write NetCDF: ", t.elapsed(), " seconds.");


/*
  // Gather back to the head node
  var GatheredArray : [Inner] real;
  GatheredArray = OutputArray;

  write_array_to_PNG(outname, GatheredArray, rgb_ptr, t);

  writeln("Elapsed time to write PNG: ", t.elapsed(), " seconds.");

  WriteOutput(GatheredArray, ImageSpace, offset);

  writeln("Elapsed time to write NetCDF: ", t.elapsed(), " seconds.");
*/
}

