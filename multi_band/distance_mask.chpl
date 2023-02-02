module distance_mask {

proc create_distance_mask(radius : real, dx : real, nx : int) {

  const D : domain(2, int) = {-nx..nx, -nx..nx};

  var dist : [D] real;

  var mask : [D] bool;

  for (i,j) in dist.domain do {
    dist[i,j] = dx * sqrt(i**2 + j**2);
  }

  // Using < here instead of <= because <= leaves only one point at the edge of the
  // domain, and it becomes difficult to define the left and right masks in a sensible way.
  mask = (dist < radius);

  var mask_size = (+ reduce mask);

  return (mask, mask_size);
}

} // distance_mask

