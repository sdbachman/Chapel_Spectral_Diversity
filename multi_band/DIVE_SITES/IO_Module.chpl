use NetCDF.C_NetCDF;
use Time;

proc CreateNetCDF(filename : string, ImageSpace : ?, out varid_out) {
  /* IDs for the netCDF file, dimensions, and variables. */
    var ncid, x_dimid, y_dimid : c_int;
    var x_varid, y_varid : c_int;
    var varid : c_int;

    var ndims : int = 2;
    var dimids: [0..#ndims] c_int;

    var shape = ImageSpace.shape;

    var yo : [1..shape[0]] real(32);
    var xo : [1..shape[1]] real(32);
    for i in 1..shape[0] {
      yo[i] = i;
    }
    for i in 1..shape[1] {
      xo[i] = i;
    }

    var yName = "y";
    var xName = "x";

  /* Create the file. */
    extern proc nc_create(path : c_string, cmode : c_int, ncidp : c_ptr(c_int)) : c_int;
    nc_create( filename.c_str(), NC_CLOBBER, c_ptrTo(ncid));

  /* Define the dimensions. The record dimension is defined to have
     unlimited length - it can grow as needed. In this example it is
     the time dimension.*/
    extern proc nc_def_dim(ncid : c_int, name : c_string, len : c_size_t, idp : c_ptr(c_int)) : c_int;
    nc_def_dim(ncid, yName.c_str(), shape[0] : c_size_t, y_dimid);
    nc_def_dim(ncid, xName.c_str(), shape[1] : c_size_t, x_dimid);

  /* Define the coordinate variables. */
    extern proc nc_def_var(ncid : c_int, name : c_string, xtype : nc_type, ndims : c_int, dimidsp : c_ptr(c_int), varidp : c_ptr(c_int)) : c_int;
      nc_def_var(ncid, yName.c_str(), NC_FLOAT, 1 : c_int, c_ptrTo(y_dimid), c_ptrTo(y_varid));
      nc_def_var(ncid, xName.c_str(), NC_FLOAT, 1 : c_int, c_ptrTo(x_dimid), c_ptrTo(x_varid));

  /* Assign units attributes to coordinate variables. */
  //  att_text = "meters";
  //  nc_put_att_text(ncid, y_varid, "units".c_str(), att_text.numBytes : c_size_t, att_text.c_str());
  //  att_text = "meters";
  //  nc_put_att_text(ncid, x_varid, "units".c_str(), att_text.numBytes : c_size_t, att_text.c_str());

  /* The dimids array is used to pass the dimids of the dimensions of
     the netCDF variables. In C, the unlimited dimension must come first on the list of dimids. */
    dimids[0] = y_dimid;
    dimids[1] = x_dimid;

  /* Define the netCDF variable. */
    nc_def_var(ncid, "spectral_diversity".c_str(), NC_FLOAT, ndims : c_int, c_ptrTo(dimids[0]), c_ptrTo(varid));
    varid_out = varid : int;

  /* Assign units attributes to the netCDF variables. */

    /*
    var maxval = max reduce(arr_out);
    var valid_range = "(0, " + (maxval : string) + ")";
    //nc_put_att_text(ncid, varid, "units".c_str(), units.numBytes : c_size_t, units.c_str());
    nc_put_att_text(ncid, varid, "Range".c_str(), valid_range.numBytes : c_size_t, valid_range.c_str());
    */

    //int nc_def_var_fill(int ncid, int varid, int no_fill, const void* fill_value)
    extern proc nc_def_var_fill(ncid : c_int, varid : c_int, no_fill : c_int, fill_value : c_ptr(c_float));
    var fv = -999 : real(32);
    nc_def_var_fill(ncid, varid, 0, c_ptrTo(fv));

  /* End define mode. */
    nc_enddef(ncid);

  /* Write the coordinate variable data. */
    extern proc nc_put_var_float(ncid : c_int, varid : c_int, op : c_ptr(c_float)) : c_int;
    nc_put_var_float(ncid, y_varid, c_ptrTo(yo[0]));
    nc_put_var_float(ncid, x_varid, c_ptrTo(xo[0]));

    nc_close(ncid);
}

var x: atomic int;

proc WriteOutput(filename : string, ref arr_out: [?D] real(64), varid_in : int, left_edge_row : [], left_edge_col : [], right_edge_col : []) {

  x.waitFor(here.id%numLocales);
  writeln("Starting on ", here.id);

  var ncid : c_int;
  var varid = varid_in : c_int;

  extern proc nc_open(path : c_string, mode : c_int, ncidp : c_ptr(c_int)) : c_int;
  nc_open( filename.c_str() , NC_WRITE, c_ptrTo(ncid));

  extern proc nc_put_var1_double(ncid: c_int, varid: c_int, indexp : c_ptr(c_size_t), op : c_ptr(c_double));
  for (i,j,k) in zip(left_edge_row, left_edge_col, right_edge_col) {
    var start = (i,j);
    var start_c : [0..#start.size] c_size_t;
    for m in 0..<start.size {
      start_c[m] = start[m] : c_size_t;
    }
    nc_put_var1_double(ncid, varid, c_ptrTo(start_c), c_ptrTo(arr_out[start]));


    for l in (j+1)..k do {
      var start = (i,l);
      var start_c : [0..#start.size] c_size_t;
      for m in 0..<start.size {
        start_c[m] = start[m] : c_size_t;
      }
      nc_put_var1_double(ncid, varid, c_ptrTo(start_c), c_ptrTo(arr_out[start]));
    }
  }
  nc_close(ncid);

  const inc = (x.read() + 1) % numLocales;
  x.write(inc);
}

inline proc tuplify(x) {
  if isTuple(x) then return x; else return (x,);
}

var y: atomic int;

proc WriteBackground(filename : string, ref arr_out: [?D] real(32), varid_in : int) {

  y.waitFor(here.id%numLocales);
  //writeln("Starting on ", here.id);

  var ncid : c_int;
  var varid = varid_in : c_int;

  extern proc nc_open(path : c_string, mode : c_int, ncidp : c_ptr(c_int)) : c_int;
  nc_open( filename.c_str() , NC_WRITE, c_ptrTo(ncid));

  writeln(arr_out.shape);

  var start = (0,0);
  var count = (D.shape[0], D.shape[1]);
  writeln(count);

  var start_c : [0..1] c_size_t;
  var count_c : [0..1] c_size_t;

  for i in 0..1 {
    start_c[i] = start[i] : c_size_t;
    count_c[i] = count[i] : c_size_t;
  }

  extern proc nc_put_vara_float(ncid : c_int, varid : c_int, startp : c_ptr(c_size_t), countp : c_ptr(c_size_t), op : c_ptr(c_float)) : c_int;
  nc_put_vara_float(ncid, varid, c_ptrTo(start_c), c_ptrTo(count_c), c_ptrTo(arr_out[start]));

  nc_close(ncid);

  const inc = (y.read() + 1) % numLocales;
  y.write(inc);
}
