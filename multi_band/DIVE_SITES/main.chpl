use distance;
//use distance_mask;
use distance_mask2;
use IO_Module;
use IO;
use BlockDist;
use Time;
use AutoMath;
use LinearAlgebra;
use IO.FormattedIO;
use AllLocalesBarriers;

/* Command line arguments. */
config const in_name : string;
config const window_size : real(32);
config const dive_site_window_size : real(32);
config const dx : real(32) = 5.0;

proc convolve_and_calculate(Array: [] real(32), dive_site : [1..2] int(16), dswr_D : ?, B : ?, Output: [] real(64), locL:?, locC:?, locR :?, Mask_Size : real(32), t: stopwatch) : [] {

  var left_edge_row = dive_site[1] + dswr_D;
  var left_edge_col = dive_site[2] - B;
  var right_edge_col = dive_site[2] + B;

  const bs = 1;
  const be = 5;

  forall (i,j,k) in zip(left_edge_row, left_edge_col, right_edge_col) {

    var tmpLL : real(64) = 0;
    var tmpLC : real(64) = 0;
    var tmpLR : real(64) = 0;
    var tmpCC : real(64) = 0;
    var tmpCR : real(64) = 0;
    var tmpRR : real(64) = 0;
    var keep  : real(64) = 0;
    var discard : real(64) = 0;
    var current : real(64) = 0;

    calc_distance(Array, locL, locL, tmpLL, bs, be, i, j);
    calc_distance(Array, locL, locC, tmpLC, bs, be, i, j);
    calc_distance(Array, locL, locR, tmpLR, bs, be, i, j);
    calc_distance(Array, locC, locC, tmpCC, bs, be, i, j);
    calc_distance(Array, locC, locR, tmpCR, bs, be, i, j);
    calc_distance(Array, locR, locR, tmpRR, bs, be, i, j);

    keep = tmpLL + tmpCC + tmpRR + 2*(tmpLC + tmpLR + tmpCR);
    Output[i,j] = keep / (Mask_Size**2);
    discard = tmpLL + 2*(tmpLC + tmpLR);

    for l in (j+1)..k do {

      tmpLL = 0;
      tmpLC = 0;
      tmpLR = 0;
      tmpCC = 0;
      tmpCR = 0;
      tmpRR = 0;

      calc_distance(Array, locL, locL, tmpLL, bs, be, i, l);
      calc_distance(Array, locL, locC, tmpLC, bs, be, i, l);
      calc_distance(Array, locL, locR, tmpLR, bs, be, i, l);
      calc_distance(Array, locC, locR, tmpCR, bs, be, i, l);
      calc_distance(Array, locR, locR, tmpRR, bs, be, i, l);

      current = tmpRR + 2*(tmpLR + tmpCR);
      keep = keep + current - discard;
      Output[i,l] = keep / (Mask_Size**2);
      //Output[i,l] = Output[i,l-1] + current - discard;
      discard = tmpLL + 2*(tmpLC + tmpLR);

      if (Output[i,l] < 0) {
        writeln("Negative output detected.");
        exit();
      }

    }
  }
}


proc main(args: [] string) {

  var t : stopwatch;
  t.start();

  ////////////////////////////////////////////////////////////////////////////////
  // Gather input variables from command line
  const window_radius = (sqrt(window_size) / 2 / dx) : int;
  writeln("Window has a radius of ", window_radius, " points.");

  const dive_site_window_radius = (sqrt(dive_site_window_size) / 2 / dx) : int;
  writeln("Dive site window has a radius of ", dive_site_window_radius, " points.");

  /////////////////////////////////////////////////////////////////////////////////
  // Open the text file that contains the number of rows and columns for the image.
  const image_infile = open(in_name + ".txt", iomode.r);
  const image_reader = image_infile.reader();
  // Read the number of rows and columns in the array in from the file.
  const rows = image_reader.read(int);
  const cols = image_reader.read(int);

  const ImageSpace = {0..<rows, 0..<cols, 1..5};

  // Close the file.
  image_reader.close();
  image_infile.close();
  ////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////
  // Open the text file that contains the dive sites
  const dive_infile = open(in_name + "_dive_sites.txt", iomode.r);
  const dive_reader = dive_infile.reader();
  // Read the number of rows and columns in the array in from the file.
  const num_dive_sites = dive_reader.read(int);
  const dive_sites_domain = {1..num_dive_sites};

  // Declare an array of the specified dimensions.
  var dive_sites: [1..num_dive_sites, 1..2] int(16);

  // Read in the array (row-major order is used for whole-array reads
  // like this).
  dive_reader.read(dive_sites);

  // Close the file.
  dive_reader.close();
  dive_infile.close();

  // Create Block distribution of dive sites
  const D = dive_sites_domain dmapped Block(dive_sites_domain);
  ////////////////////////////////////////////////////////////////////////////////

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

  ////////////////////////////////////////////////////////////////////////////////

  coforall loc in Locales do on loc {

    var t2 : stopwatch;
    t2.start();

    // Create distance mask
    const wr = (sqrt(window_size) / 2) : int;
    const wnx = (wr / dx) : int;
    const (LeftMask, CenterMask, RightMask, Mask_Size) = create_distance_mask(wr, dx, wnx);
    const locLeftMaskDomain = LeftMask.domain;
    const locCenterMaskDomain = CenterMask.domain;
    const locRightMaskDomain = RightMask.domain;

    ////////////////////////////////////////////////////////
    // Make iterator for dive site window
    const tmp = {-dive_site_window_radius..dive_site_window_radius};
    const dswr_D : [tmp] int = tmp;
    const B = round(sqrt( dive_site_window_radius**2 - dswr_D**2 )) : int;

    ////////////////////////////////////////////////////////

    const locImageSpace = ImageSpace;
    var locImage : [locImageSpace] real(32);
    var OutputArray : [0..<rows, 0..<cols] real(64);

    // Read in array
    var f = open(in_name + ".bin", iomode.r);
    var r = f.reader(kind=ionative);

    for i in locImageSpace.dim[0] {
      for j in locImageSpace.dim[1] {
        for k in locImageSpace.dim[2] {
          var tmp : real(32);
          r.readBinary(tmp);
          locImage[i,j,k] = tmp;
        }
      }
    }
    r.close();

    writeln("Elapsed time to read file on ", here.name, ": ", t.elapsed(), " seconds.");
    ////////////////////////////////////////////////////////

    for i in D.localSubdomain() {
      //writeln("Locale ", here.name, " owns ", dive_sites[i,..]);
      var dive_site = dive_sites[i,..];

      // This barrier makes sure Locale 0 prints everyone's reports before starting;
      // prevents Locale 0 from forcing other Locales to wait before they can start
      // their own convolve_and_calculate
      allLocalesBarrier.barrier();

      convolve_and_calculate(locImage, dive_site, dswr_D, B, OutputArray,locLeftMaskDomain, locCenterMaskDomain, locRightMaskDomain,Mask_Size, t2);
      t2.stop();

      writeln("Elapsed time to calculate Q on ", here.name, ": ", t2.elapsed(), " seconds.");

      var left_edge_row = dive_site[1] + dswr_D;
      var left_edge_col = dive_site[2] - B;
      var right_edge_col = dive_site[2] + B;

      allLocalesBarrier.barrier();
      WriteOutput(out_file, OutputArray, varid, left_edge_row, left_edge_col, right_edge_col);
    }
  }

  writeln("Elapsed time to write NetCDF: ", t.elapsed(), " seconds.");

}










////// Keeping this here just in case //////
/*    if (here.id == 0) {
      var yay = locImage[..,..,1];
      WriteBackground(out_file, yay, varid);
    }
*/
