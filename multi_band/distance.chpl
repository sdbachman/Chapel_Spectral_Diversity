

proc calc_distance(Array : [] real, MaskDomain1 : ?, MaskDomain2 : ?, inout Output_tmp : real, bs : int, be : int, i : int, j : int, dist_type : string) {

for (k,l) in MaskDomain1 {
  for (m,n) in MaskDomain2 {
  var dist : real = 0;

  if (dist_type == "Euclidean") {
    var tmp : real = 0;
    for p in bs..be {
      tmp += (Array[p,i+k,j+l]-Array[p,i+m,j+n])**2;
    }
    dist = sqrt(tmp);
  }
  else if (dist_type == "Manhattan") {
    var tmp : real = 0;
    for p in bs..be {
      tmp += abs(Array[p,i+k,j+l]-Array[p,i+m,j+n]);
    }
    dist = tmp;
  }
  else if (dist_type == "Canberra") {
    var tmp : real = 0;
    for p in bs..be {
      tmp += abs(Array[p,i+k,j+l]-Array[p,i+m,j+n]) / ( abs(Array[p,i+k,j+l]) + abs(Array[p,i+m,j+n]) ) ;
    }
    dist = tmp;
  }

  Output_tmp += dist;
}
}

}
