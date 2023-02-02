use IO;

var nx = 5;
var ny = 8;

var arr : [1..nx, 1..ny] real;

var f = open("/glade/work/bachman/Jupyter_Notebooks/Coral_Reef_Alliance/spectral_diversity/test.bin", iomode.r);
var r = f.reader(kind=ionative);
for i in 1..nx {
  for j in 1..ny {
    var tmp : real;
    r.readBinary(tmp);
    arr[i,j] = tmp;
  }
}
r.close();

writeln(arr);


writeln(arr.shape[0]);
