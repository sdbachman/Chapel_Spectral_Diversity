use png;
use distance_mask;
use distance;
use BlockDist;
use Time;
use IO_Module;
use IO;
use AutoMath;
use LinearAlgebra;

/* Command line arguments. */
config const in_array : string;               /* name of array to read in */
config const in_image : string;                /* name of PNG file to read */
config const outname : string;               /* name of PNG file to write at the end */
config const window_size : real;                  /* the desired area of the neighborhood (in meters^2) */
config const dx : real;                      /* the resolution of the raster image (in meters) */
config const metric : string;                /* Choose which distance metric to use */

var bs = 1;
var be = 5;

proc convolve_and_calculate(Array: [] real, centerPoints : ?, locL : ?, locC : ?, locR : ?, Output: [] real, Mask_Size : int, dist_type, t: stopwatch) : [] {

  var first_point = centerPoints.first[1];
  var last_point = centerPoints.last[1];

  forall i in centerPoints[..,first_point] {

    var tmpLL : real = 0;
    var tmpLC : real = 0;
    var tmpLR : real = 0;
    var tmpCC : real = 0;
    var tmpCR : real = 0;
    var tmpRR : real = 0;

    calc_distance(Array, locL, locL, tmpLL, bs, be, i, first_point, dist_type);
    calc_distance(Array, locL, locC, tmpLC, bs, be, i, first_point, dist_type);
    calc_distance(Array, locL, locR, tmpLR, bs, be, i, first_point, dist_type);
    calc_distance(Array, locC, locC, tmpCC, bs, be, i, first_point, dist_type);
    calc_distance(Array, locC, locR, tmpCR, bs, be, i, first_point, dist_type);
    calc_distance(Array, locR, locR, tmpRR, bs, be, i, first_point, dist_type);

    Output[i,first_point] = (tmpLL + tmpCC + tmpRR + 2*(tmpLC + tmpLR + tmpCR)) / (Mask_Size**2);
    var prev = tmpCC + tmpRR + 2*tmpCR;

    for j in (first_point+1)..last_point do {

      tmpLL = 0;
      tmpLC = 0;
      tmpLR = 0;
      tmpCR = 0;
      tmpRR = 0;

      calc_distance(Array, locL, locL, tmpLL, bs, be, i, j, dist_type);
      calc_distance(Array, locL, locC, tmpLC, bs, be, i, j, dist_type);
      calc_distance(Array, locL, locR, tmpLR, bs, be, i, j, dist_type);
      calc_distance(Array, locC, locR, tmpCR, bs, be, i, j, dist_type);
      calc_distance(Array, locR, locR, tmpRR, bs, be, i, j, dist_type);

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

  // Read in PNG
  var (rgb_ptr, Image) = load_PNG_into_array(in_image, t);
  const ImageSpace = Image.domain;
  writeln("ImageSpace is ", ImageSpace);
  writeln("Elapsed time to read into array: ", t.elapsed(), " seconds.");

  var x = Image.shape[0];
  var y = Image.shape[1];
  var Array : [1..5,1..x,1..y] real;

  // Read in array
  var f = open(in_array, iomode.r);
  var r = f.reader(kind=ionative);
  for i in 1..5 {
    for j in 1..x {
      for k in 1..y {
        var tmp : real;
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

  writeln("Elapsed time at start of coforall loop: ", t.elapsed(), " seconds.");

  writeln("Starting coforall loop.");

  coforall loc in Locales do on loc {

    // If I put "create_distance_mask" inside this loop I need to declare local copies of these variables,
    // otherwise it seems like Chapel will have to do a ton
    // of cross-locale calls to access these variables. This seems to double the amount of time it
    // takes to run through the coforall loop for all non-head locales!

    const dist_type = metric;

    const loc_Mask_Size = Mask_Size;

    const locArrayDomain = Array.domain;
    const locArray : [locArrayDomain] Array.eltType = Array;

    const locLeftMaskDomain = LeftMask.domain;
    const locCenterMaskDomain = CenterMask.domain;
    const locRightMaskDomain = RightMask.domain;

    convolve_and_calculate(locArray, D.localSubdomain(), locLeftMaskDomain, locCenterMaskDomain, locRightMaskDomain, OutputArray, loc_Mask_Size, dist_type, t);
  }


  writeln("Elapsed time to finish coforall loop: ", t.elapsed(), " seconds.");

  // Gather back to the head node
  var GatheredArray : [Inner] real;
  GatheredArray = OutputArray;

  write_array_to_PNG(outname, GatheredArray, rgb_ptr, t);

  writeln("Elapsed time to write PNG: ", t.elapsed(), " seconds.");

  WriteOutput(GatheredArray, ImageSpace, offset);

  writeln("Elapsed time to write NetCDF: ", t.elapsed(), " seconds.");
}

