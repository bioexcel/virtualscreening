# configCMIP_titration
# CMIP configuration file to run a Structure Titration.
sub configCMIP_titration{

   my ($titwat, $titip, $titim)=@_;
   $titwat=0 if (! defined $titwat);
   $titip=0 if (! defined $titip);
   $titim=0 if (! defined $titim);

   my $gr=CMIP::Grid->new();
   $gr->readGrid(2);
   $gr->perfill(0.7);
   $gr->int(0.8,0.8,0.8);

   my $inpP=CMIP::InputParams->new ("Titration", $gr);
   $inpP->addKeyword('tipcalc' =>1,
                     'titration'=>1,
                     'readgrid'=>2,
                     'inifoc'  =>2,
                     'cutfoc'  =>-0.5,
                     'focus'   =>1,
                     'ninter'  =>10,
                     'dields'  =>2,
                     'clhost'  =>1,
                     'calcgrid'=>1,
                     'irest'   =>0,
                     'orest'   =>0,
                     'coorfmt' =>2,
                     'titwat'  =>$titwat,
                     'titip'   =>$titip,
                     'titim'   =>$titim,
                      );
    return $inpP;
}


# configCMIP_interaction
# CMIP configuration file to run an Electrostatic Interaction Potential calculation.
sub configCMIP_interaction{

    my @cen=@_;

    my $gr=CMIP::Grid->new();
    $gr->readGrid(0);
    $gr->int(0.5,0.5,0.5);
    $gr->cen($cen[0],$cen[1],$cen[2]);
    $gr->dim(40,40,40);

    my $gr0=CMIP::Grid->new();
    $gr0->readGrid(2);
    $gr0->int(1.5,1.5,1.5);
    $gr0->perfill(0.6);

    my $inpP=CMIP::InputParams->new ("AA-Neutral Prot Interaction", $gr, $gr0);
    $inpP->addKeyword('tipcalc'=>3,
                      'irest'=>0,
                      'orest'=>0,
                      'cutvdw'=>8.,
                      'cutelec'=>10.,
                      'coorfmt'=>2,
                      'calcgrid'=>1,
                      'dields'=>2,
                      'outegrid'=>1,
                      'fvdw'=>0.8,
                      'pbelec'=>1
                      );
    $inpP;
}


# configCMIP_PBsolvation
# CMIP configuration file to run an Poisson-Boltzmann Solvation calculation.
sub configCMIP_PBsolvation{

    my $gr=CMIP::Grid->new();
    $gr->readGrid(2);
    $gr->perfill(0.95);
    $gr->int(0.5,0.5,0.5);

    my $gr0=CMIP::Grid->new();
    $gr0->readGrid(2);
    $gr0->perfill(0.7);
    $gr0->int(1.5,1.5,1.5);

    my $inpP=CMIP::InputParams->new ("AA-Neutral Prot PB Solvation", $gr, $gr0);
    $inpP->addKeyword('tipcalc'=>0,
                      'irest'=>0,
                      'fullrst'=>0,
                      'coorfmt'=>2,
                      'calcgrid'=>1,
                      'carmip'=>1,
                      'orest'=>0,
                      'rstonly'=>0,
                      'pbelec'=>1,
                      'pbinic'=>2,
                      'pbfocus'=>1,
                      'solvenergy'=>1
                      );
    $inpP;
}
