#! /usr/bin/perl -w

if ($#ARGV < 0) {
  print STDOUT "\nUsage: mesh2ribwireframe meshfile\n\n";
  exit(1);
}

$meshfile=$ARGV[0];

open(MESHFILE,"< ".$p_space_file);
while(defined($line=<MESHFILE>)) {
  chomp($line);
  @tokens=split(' ',$line);
  if (($tokens[0] cmp "Vertex")==0) {
    $verts[$#verts+1]
  }
}
close(MESHFILE);

$line=<PSPACE>;
chomp($line);
@barange=parsesequence($line);

$line=<PSPACE>;
chomp($line);
@Rrange=parsesequence($line);

$line=<PSPACE>;
chomp($line);
@kprange=parsesequence($line);

$line=<PSPACE>;
chomp($line);
@Drange=parsesequence($line);

close(PSPACE);


print STDOUT "[ @tfrange ]\n";
print STDOUT "[ @barange ]\n";
print STDOUT "[ @Rrange ]\n";
print STDOUT "[ @kprange ]\n";
print STDOUT "[ @Drange ]\n";

$prefix="$runtype"."_params";

@scriptname=();
#print STDOUT "\nEnter the name of the first PST script (for mcell): ";
$scriptname[0] = $runtype.".p1.".$p_space_file;
open(FILE1,"> ".$scriptname[0]);

#print STDOUT "\nEnter the name of the second PST script (for avg_dat): ";
$scriptname[1] = $runtype.".p2.".$p_space_file;
open(FILE2,"> ".$scriptname[1]);

#print STDOUT "\nEnter the name of the third PST script (for blip): ";
$scriptname[2] = $runtype.".p3.".$p_space_file;
open(FILE3,"> ".$scriptname[2]);

$count=0;
$totaltaskload=0;
@taskloadarray=();
@outputnameprefix=();
for ($i1=0;$i1<=$#Rrange;$i1++) {
  for ($i2=0;$i2<=$#kprange;$i2++) {
    for ($i3=0;$i3<=$#Drange;$i3++) {
      for ($i4=0;$i4<=$#barange;$i4++) {
	for ($i5=0;$i5<=$#tfrange;$i5++) {
          $paramfilename="$prefix.".$count.".mcell";
          print STDOUT "Generating $paramfilename\n";
          $taskload=calcParams($paramfilename,$runtype,$tfrange[$i5],$barange[$i4],$Rrange[$i1],$kprange[$i2],$Drange[$i3]);
          $taskload = sprintf "%d",$taskload/1000;
          $totaltaskload+=$taskload;
          $taskloadarray[$#taskloadarray+1]=$taskload;
  
	  $outputnameprefix[$count]=get_output_prefix(
             $runtype,$tfrange[$i5],$barange[$i4],$Rrange[$i1],$kprange[$i2],
             $Drange[$i3]);

          for ($seed=0;$seed<=$#seeds;$seed++) {
            $taskid=($count)*($#seeds+1)+$seed;
            print FILE1 "$taskid $taskload mcell $prefix.$count.mcell $seeds[$seed]\n";
          }

          print FILE2 "$count 1 avg_dat ";
          for ($seed=0;$seed<=$#seeds;$seed++) {
            print FILE2 "$outputnameprefix[$count]"."seed_$seeds[$seed].A2Ro ";
          }
          print FILE2 "$outputnameprefix[$count]"."avg.A2Ro\n";

          print FILE3 "$count 1 blip $outputnameprefix[$count]"."avg.A2Ro";
          print FILE3 " > $outputnameprefix[$count]"."blip.A2Ro\n";

	  $count++;
        }
      }
    }
  }
}
close(FILE1);
close(FILE2);
close(FILE3);
$numberoftasks=$taskid+1;

@taskloadarray = sort {$a <=> $b} @taskloadarray;

#$loadaccum=0;
#for ($i=0;$i<=$#taskloadarray;$i++) {
#  $loadaccum+=$#seeds*$taskloadarray[$i];
#  if ($loadaccum<=$totaltaskload/2) {
#    $taskhalflife=$taskloadarray[$i];
#  }
#}

$totaltaskload=($#seeds+1)*$totaltaskload;
print STDOUT "total task load = $totaltaskload\n";
print STDOUT "total number of tasks = $numberoftasks\n";
printf STDOUT "avg load per task = %.9g\n",$totaltaskload/$numberoftasks;
#printf STDOUT "task 1/2 life = %.9g\n",$taskhalflife;

open(NODELIST,"/usr/local/bin/LOADL_NODE_LIST.pl |");
#open(NODELIST,"< ".$nodelistfile);
@nodelist=();
while(defined($line=<NODELIST>)) {
  chomp($line);
  @nodelist=split(/ +/,$line);
}
close(NODELIST);

$numberofnodes=$#nodelist+1;
$numberofcpus=$numberofnodes*$cpuspernode;
print STDOUT "number of nodes = $numberofnodes\n";
$loadpernode=$totaltaskload/$numberofnodes;
$loadpercpu=$totaltaskload/$numberofcpus;
print STDOUT "load per node = $loadpernode\n";
print STDOUT "load per cpu = $loadpercpu\n";

@cpunamelist=();
@cpuloadlistsort=();
keys(%cpuloading)=1000;
for ($nodenum=0;$nodenum<$numberofnodes;$nodenum++) {
  for ($cpunum=0;$cpunum<$cpuspernode;$cpunum++) {
    $cpuid=($cpuspernode*$nodenum)+$cpunum;
    $cpunamelist[$cpuid]=sprintf "$nodelist[$nodenum] $cpunum"; 
    $cpuloadlistsort[$cpuid]=$cpuid;
    $cpuloading{$cpuid}=0;
  }
}

@tasklist=();
keys(%taskloading)=1000;
open(FILE1,"< ".$scriptname[$scriptnum]);
while(defined($line=<FILE1>)) {
  chomp($line);
  @taskline=split(/ +/,$line);
  $tasklist[$#tasklist+1]=$line;
  $taskloading{$taskline[0]}=$taskline[1];
}
close(FILE1);
$numberoftasks=$#tasklist+1;
@taskloadlistsort = sort { $taskloading{$b} <=> $taskloading{$a} } keys %taskloading;

$loadaccum=0;
$cpuindex=0;
keys(%cputasklist)=1000;
@taskqueue=();
for ($i=0;$i<$numberoftasks;$i++) {
  @taskline=split(/ +/,$tasklist[$taskloadlistsort[$i]]);
  $taskload=$taskline[1];
  $loadaccum+=$taskload;
  $cpuid=$cpuloadlistsort[$cpuindex];
  $cputasklist{$i}=$cpunamelist[$cpuid];
  $taskqueue[$i]=sprintf "$cpunamelist[$cpuid]";
  $cpuloading{$cpuid}+=$taskload;
  for ($j=0;$j<=$#taskline;$j++) {
    $taskqueue[$i]=sprintf "%s %s",$taskqueue[$i],$taskline[$j];
  }
#  if ($loadaccum>$loadpercpu) {
#    $loadaccum=0;
    $cpuindex++;
    if ($cpuindex == $numberofcpus) {
      $cpuindex=0;
      @cpuloadlistsort = sort { $cpuloading{$a} <=> $cpuloading{$b} or $a <=> $b; } keys %cpuloading;
    }
#  }
}

@cputasklistsort = sort { $cputasklist{$a} cmp $cputasklist{$b} } keys %cputasklist;
open(FILE2,"> ".$scriptname[$scriptnum].".tasklist");
for ($i=0;$i<$numberoftasks;$i++) {
  print FILE2 "$taskqueue[$cputasklistsort[$i]]\n";
}
close(FILE2);

print STDOUT "*****\n";
foreach $key (sort { $cpuloading{$a} <=> $cpuloading{$b} or $a <=> $b; } keys %cpuloading) {
  printf STDOUT "%4d %s\n", $key, $cpuloading{$key};
}
print STDOUT "*****\n";

if ($starttask == 1) {
  $taskfile=$scriptname[$scriptnum].".tasklist";
  print STDOUT "executing /paci/salk/u12045/scripts/start_task2 $taskfile\n";
  `/paci/salk/u12045/scripts/start_task2 $taskfile`;
#  print STDOUT "executing /paci/salk/u12045/scripts/start_task $taskfile\n";
#  `/paci/salk/u12045/scripts/start_task $taskfile`;
}

exit(0);


##############
# get_output_prefix()
#
sub get_output_prefix {
  my($runtype,$tf,$ba,$R,$kp,$D)=@_;
  my($result);

# Form the output prefix;
  $result="";
  $result=$runtype."_".
          "tf_".sprintf("%.3g",$tf*1e3)."_".
          "ba_".sprintf("%.3g",$ba)."_".
          "R_".sprintf("%.3g",$R)."_".
          "kp_".sprintf("%.3g",$kp*1e-7)."_".
          "D_".sprintf("%.3g",$D*1e6)."_";
 
  return $result;
}


##############
# parsesequence
#
# TO DO :
#    - ERROR CHECKING !!!
#
sub parsesequence {
  my($string)=@_;
  my(@list,@tmp);

  @tmp = split(/[ ,]+/,$string);

  $count=0;
  for ($i=0;$i<=$#tmp;$i++) {
    @tokens= split(/:/,$tmp[$i]);
    $begin=$tokens[0];
    $step=$tokens[1];
    $end=$tokens[2];
    if (!defined($step)) {
      $step = 1;
    }
    if (!defined($end)) {
      $end = $begin;
    }
    chomp($begin);
    chomp($end);
    chomp($step);
    for ($j=$begin;$j<$end+$step;$j+=$step) {
      $list[$count++] = $j;
    }
  }
  return @list;
}


##############
# calcParams()
#
sub calcParams {
  my($filename,$runtype,$tf,$ba,$R,$kp,$D)=@_;
  my($alpha,$beta,$km,$achr_d,$grid_density,$total_sim_time,$PI,$Na);
  my($p_plus_target,$dt_factor,$dt_kp,$dt,$it,$iterations,$f);

  open(FILE,"> ".$filename);
  $alpha = (1+$R)/$tf;
  $beta = $ba * $alpha;
  $km = $beta/(2*$R);
  $achr_d = 7250;

  $grid_density = 10000;
  $total_sim_time=40e-3;

  $PI=3.14159265358979323846;
  $Na=6.02205e23;
  $p_plus_target=0.15;
  $dt_factor=(($p_plus_target*$Na*1e-11/$grid_density)**2)*$D/$PI;
  $dt_kp = $dt_factor/($kp**2);
  $dt=sprintf "%.2g",$dt_kp;
  
  if ($dt>=1e-6) {
    $dt=1e-6;
  }

  $it=$total_sim_time/$dt;
  $iterations=sprintf "%d",$it;

  $f=sprintf "%.9g",2e-7/$dt;
  if ($f >= 1) {
    $f=sprintf "%d",2e-7/$dt;
  }
  else {
    $f=1;
  }
  printf FILE "tf = %.9g\n",$tf;
  printf FILE "ba = %.9g\n",$ba;
  printf FILE "R = %.9g\n",$R;
  printf FILE "kp = %.9g\n",$kp;
  printf FILE "D = %.9g\n",$D;

  print FILE "run_type = \"$runtype\"\n";
  if ($runtype eq "dfp") {
    print FILE "btx_factor = 1.0\n";
    print FILE "ache_d_factor = 0.03\n";
  }
  elsif ($runtype eq "hbtx") {
    print FILE "btx_factor = 0.295\n";
    print FILE "ache_d_factor = 0.03\n";
  }
  elsif ($runtype eq "norm") {
    print FILE "btx_factor = 1.0\n";
    print FILE "ache_d_factor = 1.0\n";
  }

  printf FILE "dt = %.9g\n",$dt;
  printf FILE "ITERATIONS = %d\n",$it;
  printf FILE "TIME_STEP = %.9g\n",$dt;
  printf FILE "freq = %d*dt\n",$f;
  print FILE "INCLUDE_FILE = \"main.mdl\"\n";
  close(FILE);
 
  return $iterations;

}
