
proc calc_distance(Array : [] real(32), MaskDomain1 : ?, MaskDomain2 : ?, inout Output_tmp : real(64), bs : int, be : int, i : int, j : int) {

for (k,l) in MaskDomain1 {
  for (m,n) in MaskDomain2 {
    var dist : real(64) = 0;
    var tmp : real(64) = 0;

    for p in bs..be {
      tmp += (Array[i+k,j+l,p]-Array[i+m,j+n,p])**2;
    }
    dist = sqrt(tmp);
    Output_tmp += dist;

  }
}

}
