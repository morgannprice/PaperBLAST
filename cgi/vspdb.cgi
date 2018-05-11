#!/usr/bin/perl -w
use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use IO::Handle; # for autoflush

# CGI filtering arguments, all optional:
# org -- i.e. "Homo sapiens"
# complete -- 1 for only show proteins with no information whatsoever.
# soluble -- "yes", "partly", "no", or "" for any
# PFam -- 1 to only show proteins with a PFam domain
# ordered -- 1 to exclude disordered proteins
# n -- number of lines (default is 50)

my $filterOrg = param('org') || "";
my $filterComplete = param('complete') || 0;
my $filterSoluble = param('soluble') || "";
my $filterPFam = param('PFam') || 0;
my $filterOrder = param('ordered') || 0;
my $disorderThreshold = 0.4;
my $nToShow = param('n') || 50;

my @orgChoices = sort ("Caenorhabditis elegans", "Saccharomyces cerevisiae", "Schizosaccharomyces pombe", "Pseudomonas aeruginosa", "Chlamydia trachomatis", "Mycobacterium tuberculosis", "Drosophila melanogaster", "Arabidopsis thaliana", "Rattus norvegicus", "Escherichia coli", "Mus musculus", "Homo sapiens");

my %orgLabels = map { $_ => $_ } @orgChoices;
$orgLabels{""} = "Any";
my @orgChoices2 = @orgChoices; unshift @orgChoices2, "";
my $labelYesNo = {0=>"No",1=>"Yes"};

my $orgSelector = 
my $completeSelector = 
my $pfamSelector = 
my $disorderSelector = 
my $solubleSElector = 
my $nSelector = 

my $title = "Papers vs. PDB";
print
  header(-charset => 'utf-8'),
  start_html(-head => Link({-rel => "shortcut icon", -href => "../static/favicon.ico"}),
             -title => $title),
  h2("Papers vs. PDB: Well-studied proteins that lack structural information."),
  p(a({-href => "vspdb.cgi"}, "Papers vs. PDB"),
    "shows proteins that are linked to at least 10 papers yet they lack homologs in the",
    a({-href => "http://www.rcsb.org"}, "RSCB Protein Data Bank") . "."),
  h3("Filtering:"),
  start_form( -name => 'input', -method => 'GET', -action => 'vspdb.cgi'),
  p("Organism:", popup_menu('org', \@orgChoices2, $filterOrg, \%orgLabels)),
  p("Entire protein lacks homologs in PDB:",
    popup_menu('complete', [0,1], $filterComplete, $labelYesNo),
    "Has a hit in PFam:",
    popup_menu('PFam', [0,1], $filterPFam, {1=>"Yes", 0=>"Either"}),
    "Exclude Disordered:",
    popup_menu('ordered', [0,1], $filterOrder, $labelYesNo) ),
  p("Solubility:",
    popup_menu('soluble', ["","yes","partly","no"], $filterSoluble,
               { "" => "Any", "yes" => "Entirely soluble", "partly" => "Partly soluble", "no" => "Integral membrane protein"}) ),
  p("Proteins to show:", popup_menu('n', [20,50,100,200,500], $nToShow)),
  p(submit('Go'), reset()),
  end_form();

my $tsv = "../vspdb/unhit.tsv";
open(TSV, "<", $tsv) || die "Cannot read $tsv\n";
my $header = <TSV>;
chomp $header;
my @cols = split /\t/, $header;
my %colindex = map { $cols[$_] => $_ } (0..(scalar(@cols)-1));
foreach my $field (qw{uniqId begin end organism source sequenceId desc nPapers soluble nTMH nSoluble PFams nPFam IUPredL}) {
  die "Field $field is missing from $tsv\n"
    unless exists $colindex{$field};
}

autoflush STDOUT 1; # show preliminary results
print CGI::start_table({ -cellpadding => 3, -cellspacing => 0 }), "\n";
print Tr({-valign => "top"},
         th({ -colspan => 3}, "Protein that Lacks Structural Information"),
         th(a({ -title => "The amino acid coordinates of the region of the protein that lacks homologs"}, "Region")),
         th("#Papers")) . "\n";

my %sourceToURL = ( "SwissProt" => "http://www.uniprot.org/uniprot/",
                    "SwissProt/TReMBL" => "http://www.uniprot.org/uniprot/",
                    "BRENDA" => "http://www.uniprot.org/uniprot/",
                    "MicrobesOnline" => "http://www.microbesonline.org/cgi-bin/fetchLocus.cgi?locus=",
                    "RefSeq" => "http://www.ncbi.nlm.nih.gov/protein/",
                    "metacyc" => "https://metacyc.org/gene?orgid=META&id=",
                    "ecocyc" => "https://ecocyc.org/gene?orgid=ECOLI&id=",
                    "CharProtDB" => "",
                    "CAZy" => "http://www.cazy.org/search?page=recherche&lang=en&tag=4&recherche="
                  );

my $nShown = 0;
while(my $line = <TSV>) {
  chomp $line;
  my @values = split /\t/, $line;
  die "Invalid line $line\n" unless scalar(@values) == scalar(@cols);
  my %values = map { $cols[$_] => $values[$_] } (0..(scalar(@cols)-1));
  next unless $filterOrg eq "" || $values{organism} =~ m/^$filterOrg/;
  next if $filterComplete && ($values{begin} != 1 || $values{end} != $values{length});
  if ($filterSoluble eq "yes") {
    next unless $values{nTMH} == 0;
  } elsif ($filterSoluble eq "partly") {
    next unless $values{soluble} == 1;
  } elsif ($filterSoluble eq "no") {
    next unless $values{nTMH} > 1;
  }
  next if $filterPFam && $values{nPFam} == 0;
  next if $filterOrder && $values{IUPredL} > $disorderThreshold;

  my $disorder = int($values{IUPredL} * 100);

  my @pfamlinks = ();
  my @pfamextra = ();
  my %pfamShown = ();
  foreach my $pfamspec (split /,/, $values{PFams}) {
    my ($pfam, $name, $dombeg, $domend) = split /:/, $pfamspec;
    $pfam =~ s/[.]\d+$//; # from model name to pfam id, i.e. PF12820.6 to PF12820
    next if exists $pfamShown{$pfam};
    $pfamShown{$pfam} = 1;
    if (@pfamlinks >= 3) {
      push @pfamextra, "$name";
    } else {
      push @pfamlinks, a({ -href => "http://pfam.xfam.org/family/$pfam", -title => "PFam $pfam $dombeg:$domend", }, $name);
    }
  }
  push @pfamlinks, a({ -title => join(", ", @pfamextra) }, "...")
    if @pfamextra;
  my $pfamShow = @pfamlinks == 0 ? "no PFam hits" : "PFams: " . join(", ", @pfamlinks);


  my $idShort = $values{sequenceId}; $idShort =~ s/^.*:://;
  my $newline = "%0A";
  my $url_TMHMM = "http://www.cbs.dtu.dk/cgi-bin/webface2.fcgi?"
    . join("&",
           "configfile=/usr/opt/www/pub/CBS/services/TMHMM-2.0/TMHMM2.cf",
           "outform=-noshort",
           "SEQ=>${idShort}$newline$values{subseq}");
  my $showSoluble;
  if ($values{nTMH} == 0) {
    $showSoluble = a({-title => "no transmembrane helices", -href => $url_TMHMM}, "no TMH");
  } elsif ($values{soluble} == 0) {
    $showSoluble = a({-href => $url_TMHMM}, "$values{nTMH} TMH");
  } else {
    $showSoluble = a({ -title => "$values{nSoluble} soluble amino acids and $values{nTMH} transmembrane helices",
                       -href => $url_TMHMM},
                     "partly");
  }

  my $idShort2 = $idShort;
  $idShort2 =~ s/VIMSS// if $values{source} eq "MicrobesOnline";
  my $URL = "";
  if ($values{source} eq "reanno") {
    my ($orgId, $locusId) = split /:/, $idShort;
    $URL = "http://fit.genomics.lbl.gov/cgi-bin/singleFit.cgi?orgId=$orgId&locusId=$locusId";
  } elsif ($values{source} eq "REBASE") {
    $URL = "http://rebase.neb.com/rebase/enz/${idShort}.html";
  } else {
    die "No URL for $values{source}" unless exists $sourceToURL{ $values{source} };
    $URL = $sourceToURL{ $values{source} } . $idShort2
      if $sourceToURL{ $values{source} };
  }
  my $region = $values{begin} == 1 && $values{end} == $values{length} ?
    "All $values{length}" : "$values{begin}:$values{end}/$values{length}";

  my @row1 = ( td({ -colspan => 3 },
                  a({ -title => "description from $values{source}", -href => $URL }, $values{desc}) ),
               td($region),
               td({ -align => "right" },
                  a({ -title => "Number of papers", -href => "litSearch.cgi?query=" . $values{subseq} },
                    $values{nPapers} ))
             );

  my @row2 = ( td("&nbsp;"),
               td($values{organism}),
               td(small($pfamShow)),
               td(small(a({-title => "%disordered from IUPredL",
                          -href => "http://iupred.enzim.hu/pred.php?WS=500&type=long&output=graph&seq=$values{subseq}" }, "$disorder% dis."))),
               td(small($showSoluble))
             );
  my $bgcolor = $nShown % 2 ? "white" : "lightgrey";
  print
    Tr({ -valign => "top", -bgcolor => $bgcolor },@row1),
      Tr({ -valign => "top", -bgcolor => $bgcolor }, @row2),
        "\n";
  $nShown++;
  last if $nShown >= $nToShow;
}
print CGI::end_table() . "\n";
close(TSV) || die "Error reading $tsv\n";
print p(small("Only the $nShown most-cited proteins were shown.")) if $nShown >= $nToShow;

my $date = `date -r $tsv "+%B %e, %Y"`; $date =~ s/\s+$//;
my $disorderThresholdPercent = $disorderThreshold * 100;
print h3("Technical Details"),
  start_ul(),
  li("This analysis is from $date."),
  li("A protein is considered to lack structural information if at least half of the protein,",
    "constituting at least 50 amino acids, has no homolog in PDB with &ge;30% average amino acid identity."),
  li("Transmembrane helices were identified with",
     a({-href => "http://www.cbs.dtu.dk/services/TMHMM"}, "TMHMM.")),
  li("A protein is considered partly soluble if it has one or more TMHs but a region of at least 50 amino acids has none."),
  li("Average disorder was predicted using", a({-href => "http://iupred.enzim.hu/"}, "IUPredL.")),
  li("A protein is considered disordered if the averarage disorder was above ${disorderThresholdPercent}%."),
  end_ul(),
  h3("Downloads"),
  p(a({ -href => $tsv }, "The complete table"), "(tab-delimited)."),
  end_html;
