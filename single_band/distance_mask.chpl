module distance_mask {

proc create_distance_mask(radius : real, dx : real, nx : int) {

  // var nx : int = (radius / dx) : int;

  const D : domain(2, int) = {-nx..nx, -nx..nx};
  var D_left : sparse subdomain(D);
  var D_right : sparse subdomain(D);

  var dist : [D] real;

  var center_mask : [D] bool;
  var left_mask : [D_left] bool;
  var right_mask : [D_right] bool;

  // Set default value of sparse arrays to be true.
  left_mask.IRV = true;
  right_mask.IRV = true;

  for (i,j) in dist.domain do {
    dist[i,j] = dx * sqrt(i**2 + j**2);
  }

  // Using < here instead of <= because <= leaves only one point at the edge of the
  // domain, and it becomes difficult to define the left and right masks in a sensible way.
  center_mask = (dist < radius);

  // Add some sparse indices to the sparse domains
  // "left" will represent the left edge of the circle,
  // and "right" will represent the right edge

  // left
  for (i,j) in dist.domain do {
    if (center_mask[i,j] == true && (j == -nx || center_mask[i,j-1] == false)) {
      D_left += (i,j);
    }
  }

  // right
  for (i,j) in dist.domain do {
    if (center_mask[i,j] == true && (j == nx || center_mask[i,j+1] == 0)) {
      D_right += (i,j);
    }
  }

  // Remove points from the center mask that overlap with the left and right masks
  for (i,j) in D_left do {
    center_mask[i,j] = false;
  }
  for (i,j) in D_right do {
    center_mask[i,j] = false;
  }

  var mask_size = (+ reduce left_mask) + (+ reduce center_mask) + (+ reduce right_mask);

  return (left_mask, center_mask, right_mask, mask_size);
}

} // distance_mask

