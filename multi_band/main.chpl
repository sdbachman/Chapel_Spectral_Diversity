use png;
use distance_mask;
use BlockDist;
use Time;
use IO;
use AutoMath;
use LinearAlgebra;

/* Command line arguments. */
config const in_array : string;               /* name of array to read in */
config const in_image : string;                /* name of PNG file to read */
config const outname : string;               /* name of PNG file to write at the end */
config const window_size : real;                  /* the desired area of the neighborhood (in meters^2) */
config const dx : real;                      /* the resolution of the raster image (in meters) */

proc convolve_and_calculate(Array: [] real, centerPoints : ?, Mask : [] bool, MaskDomain : ?, Output: [] real, Mask_Size : int,  t: stopwatch) : [] {

  //var centerPoints2 = [ (95,123), ];

  forall (i,j) in centerPoints {
    for (k,l) in MaskDomain {
      var dist : real = 0;
      if Mask[k,l] {
        //writeln( Array[..,i,j] );
        //writeln( Array[..,i+k,j+l] );
        var tmp : real = 0;
        for m in 1..5 {
          tmp += (Array[m,i+k,j+l]-Array[m,i,j])**2;
        }
        dist = sqrt(tmp);
        //writeln(dist);
        //writeln();
      }
      Output[i,j] += dist;
    }
    Output[i,j] = Output[i,j] / Mask_Size;
    //writeln(Output[i,j]);
    //if (Output[i,j] > 250) { writeln(i : string, " ", j : string, " ", Output[i,j] : string); }
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
  var (Mask, Mask_Size) = create_distance_mask(radius, dx, nx);

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

    const loc_Mask_Size = Mask_Size;

    const locArrayDomain = Array.domain;
    const locArray : [locArrayDomain] Array.eltType = Array;

    const locMaskDomain = Mask.domain;
    const locMask : [locMaskDomain] Mask.eltType = Mask;

    convolve_and_calculate(locArray, D.localSubdomain(), locMask, locMaskDomain, OutputArray, loc_Mask_Size, t);
  }


  writeln("Elapsed time to finish coforall loop: ", t.elapsed(), " seconds.");

  // Gather back to the head node
  var GatheredArray : [Inner] real;
  GatheredArray = OutputArray;

  write_array_to_PNG(outname, GatheredArray, rgb_ptr, t);

  writeln("Elapsed time to write PNG: ", t.elapsed(), " seconds.");
}

