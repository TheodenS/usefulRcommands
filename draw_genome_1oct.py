from Bio.Seq import Seq
import re
from Bio import SeqIO
#from Bio.Alphabet import IUPACUnambiguousDNA
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
#from commahandler import tabinout
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
import GDUtilities
from Bio.Graphics.GenomeDiagram import FeatureSet as GDFeatureSet, Diagram as GDDiagram
#from Bio.Graphics.GenomeDiagram import Graph

import Bio
import sys

print sys.path
print Bio.__version__
#gbfile="/Users/security/Dropbox/isproject_results/isproject_gbfiles_all/NoAz/NC_014248.gbk"
gbfile="/Users/security/science/frankiaproj/frankiajul2014/collected_genomes_july_2014/FraDas/NC_015656.gbk"
isgbfile="/Users/security/science/frankiaproj/frankiajul2014/outputestgb/RESULTS_BADREMOVED_over500_e-6_e-5_7sep_gbfiles/FraDas/NC_015656.gbk"
#gbfile="/Users/security/science/frankiaproj/frankiajul2014/gbresults/gb/RESULTS_BADREMOVED_over500_e-6_e-5_7sep_gbfiles/FraDas/NC_015656.gbk"
gbfile_handle=open(gbfile,"r")
parsed_gb=list(SeqIO.parse(gbfile_handle,"genbank"))[0]

isgbfile_handle=open(isgbfile,"r")
isparsed_gb=list(SeqIO.parse(isgbfile_handle,"genbank"))[0]

iselementslist= ['REMBX_TRANS_SCORE1000', 'BLASTN_CLOSEIS', 'ARB_TRANS', 'GOOD_IS_AA', 'HIGH_TRANS','OVER20_GOOD_REMX_TRANS', 'REMBX_TRANS_RATIO15', 'OVER20_GOOD_LOCX_TRANS', 'OVER15_TRANS', 'ISFINDER_CLOSE', 'GOOD_NX_LOCAL_IS', 'ISFINDER', 'ISFINDER_NOFULL','BLASTN_CLOSEIS_NOHIGH']

gd_diagram = GenomeDiagram.Diagram("stuff",tracklines=False)
#gd_diagram = GenomeDiagram.Diagram("stuff",tracklines=True)
test=GenomeDiagram.GraphSet
#raw_input("wait")
#t=test.new_graph()
#test=GenomeDiagram.Graph()
#print dir(gd_diagram)
#gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features",scale_ticks=False,height=0.2)

gd_track_for_is= gd_diagram.new_track(1,name="nifs",scale_ticks=False,scale=False,height=1,graytrack=1)
track_for_track=gd_diagram.new_track(2,name="tracks",scale_ticks=True,height=2)
gd_track_for_nif = gd_diagram.new_track(3,name="nifs",scale_ticks=False,scale=False,height=1,graytrack=1)
#gd_track_for_rna = gd_diagram.new_track(3,name="rna",scale_ticks=False,scale=False,height=1)
track_for_pseudos=gd_diagram.new_track(4,name="pseudos",scale_ticks=False,height=1,scale=False)
gd_track_for_features = gd_diagram.new_track(5,name="is",scale_ticks=False,scale=False,height=1)
gd_track_for_ruler = gd_diagram.new_track(6,name="ruler",scale_ticks=True,scale=True,height=1)



def getWindowScore(totalstring,position,stretch):

    wholestring=[]

    if position+stretch>len(totalstring):

        plusposition1=totalstring[position:len(totalstring)]
        plusposition2=totalstring[0:position+stretch-len(totalstring)]
        minusposition=totalstring[position-stretch-1:position-1]
        #wholestring=minusposition+['|']+[totalstring[position-1]]+['|']+plusposition1+plusposition2
        wholestring=minusposition+[totalstring[position-1]]+plusposition1+plusposition2

    elif position-stretch<1:

        minusposition1=totalstring[0:position-1]
        minusposition2=totalstring[position+len(totalstring)-stretch-1:len(totalstring)]
        plusposition=totalstring[position:position+stretch]
        #wholestring=minusposition2+minusposition1+['|']+[totalstring[position-1]]+['|']+plusposition
        wholestring=minusposition2+minusposition1+[totalstring[position-1]]+plusposition
    else:   

        plusposition=totalstring[position:position+stretch]
        minusposition=totalstring[position-stretch-1:position-1]
        #wholestring=minusposition+['|']+[totalstring[position-1]]+['|']+plusposition
        wholestring=minusposition+[totalstring[position-1]]+plusposition

    return wholestring.count("x")/float(len(wholestring))



#gd_track2 = gd_diagram.new_track(1, greytrack=1, name="GC content")

#gd_track_for_features.scale_ticks=True

#gd_track_for_features.name="NOAZ" 

#gd_track_for_features.height=1 
#is a float denoting the height of the track, relative to other tracks on the diagram (default=1).

#gd_track_for_features.hide=0 
#is a Boolean specifying whether the track should be drawn or not (default=0). Only tracks with hide=0 are listed by the GDDiagram method get_drawn_levels().

#gd_track_for_features.greytrack=0 
#is a Boolean specifying whether the track should include a grey background (useful for delineating many closely-spaced tracks), and a set of foreground labels (default=0).
#gd_track_for_features.greytrack_labels=5 
#is an integer specifying the number of foreground labels that should be included on the track (default=5).

#gd_track_for_features.greytrack_fontsize=8
#is an integer specifying the size of font to be used on the foreground labels (default=8).

#gd_track_for_features.greytrack_font='Helvetica'
# is a string specifying the name of the font to be used for the foreground labels. Not all machines will provide the same font selection (default=Helvetica).
#gd_track_for_features.greytrack_font_rotation=0
#is an integer specifying the angle in degrees through which to rotate the foreground labels, which are, by default, radial in circular diagrams and collinear with the track for linear diagrams (deafult=0).

#gd_track_for_features.greytrack_font_colour=colors.Color(0.6,0.6,0.6)
#is a ReportLab colors.Color object defining the colour of the foreground labels (default=colors.Color(0.6,0.6,0.6)).
#gd_track_for_features.scale=1 
#is a Boolean defining whether the track will carry a scale (default=1).

#gd_track_for_features.scale_color=colors.black
#is a ReportLab colors.Color object defining the colour of the gd_track_for_features.scale (default=colors.black).

#gd_track_for_features.scale_smallticks=0.3 
#is a float describing the height of the small tick set as a
#proportion of half the track height. Positive values run upwards (linear)
#or away from the centre of the diagram (circular), while negative values
#run in the opposite direction (default=0.3).

gd_track_for_ruler.scale_largetick_interval=1000000 
#track_for_track.scale_largetick_interval=1000000 
#is an integer specifying the interval between large
#ticks as a number of bases (default=1000000).

gd_track_for_ruler.scale_smalltick_interval=100000 
track_for_track.scale_smalltick_interval=100000 
#track_for_track.scale_smalltick_interval=100000 
#is an integer specifying the interval between smallticks
#as a number of bases (default=10000).

gd_track_for_features.scale_largetick_labels=0
#is a Boolean describing whether labels marking tick
#position will be placed over every large tick (default=1).  
#gd_track_for_features.scale_smalltick_labels=0 
gd_track_for_ruler.scale_smalltick_labels=10 
gd_track_for_ruler.scale_tick_labels=0
#is a Boolean describing whether labels marking tick
#position will be places over every small tick (default=0).




    #countposition=0
    #number_of_points=len(parsed_gb)/point_every

    #dataset=[]

    #org_genome_path=origin_gb_file_dir+orgabb+"/"+largest_contig
    ##raw_input(org_genome_path)
    #totstring=getTotalString(origin_gb_file_dir+orgabb,largest_contig,"CDS")

    #while countposition<len(parsed_gb):
    #    dataset.append((countposition,getWindowScore(totstring,countposition,windowsize)))
    #    countposition+=point_every



pseudo_set=track_for_pseudos.new_set()
gd_feature_set = gd_track_for_features.new_set()
#rna_set = gd_track_for_rna.new_set()
gd_feature_set2 = track_for_track.new_set("graph")
nif_feature_set = gd_track_for_nif.new_set()
is_feature_set = gd_track_for_is.new_set()

print "-------------------"
#print dir(gd_feature_set)
print "-------------------"
#raw_input("wait")
#calc_at_skew
graphdata = GDUtilities.gc_content(parsed_gb.seq, 2000)

totalstring=["-"]*len(parsed_gb)

for isgbfeature in isparsed_gb.features:
    is_feature_set.add_feature(isgbfeature, name="test",color=colors.HexColor('#FF6699'),label_angle=0,label_position="middle",label=False,sigil="ARROW",arrowshaft_height=1)

for gbfeature in parsed_gb.features:
    # Hypothetical proteins, gray  
    if gbfeature.type=="CDS" and gbfeature.qualifiers["product"][0]=="hypothetical protein":
        pseudo_set.add_feature(gbfeature,color=colors.HexColor('#808080'),name="test2",label_angle=0,label_position="middle",label=False,sigil="ARROW",arrowshaft_height=1)
    
    mark_color=colors.HexColor('#33CCFF')
   
    #identified proteins
    if gbfeature.type=="CDS" and gbfeature.qualifiers["product"][0]!="hypothetical protein":
        gd_feature_set.add_feature(gbfeature, name="test",color=mark_color,label_angle=0,label_position="middle",label=False,sigil="ARROW",arrowshaft_height=1)

    mark_color=colors.red
    #nitrogenase proteins
    if gbfeature.type=="CDS" and "nitrogenase" in gbfeature.qualifiers["product"][0].lower():
        nif_feature_set.add_feature(gbfeature, name="test",color=colors.HexColor('#5566FF'),label_angle=0,label_position="middle",label=False,sigil="ARROW",arrowshaft_height=1)

    if gbfeature.type=="CDS" and "rieske" in gbfeature.qualifiers["product"][0].lower():
        nif_feature_set.add_feature(gbfeature, name="test",color=colors.HexColor('#CC66FF'),label_angle=0,label_position="middle",label=False,sigil="ARROW",arrowshaft_height=1)


    mark_color=colors.HexColor('#33CC33')
    # RNA
    #if gbfeature.type=="tRNA" or gbfeature.type=="rRNA" or gbfeature.type=="ncRNA" or gbfeature.type=="tmRNA":
    #    rna_set.add_feature(gbfeature, name="test",color=mark_color,label_angle=0,label_position="middle",label=False,sigil="ARROW",arrowshaft_height=1)

    # count coding 
    if gbfeature.type=="CDS" or gbfeature.type=="tRNA" or gbfeature.type=="rRNA" or gbfeature.type=="ncRNA" or gbfeature.type=="tmRNA":
        temp_locstart=int(gbfeature.location.start)+1
        temp_locend=int(gbfeature.location.end)
        locstart=min(temp_locstart,temp_locend)
        locend=max(temp_locstart,temp_locend)
        totalstring[locstart-1:locend]=["x"]*(1+locend-locstart)

graph = gd_feature_set2.new_graph(graphdata, 'GC content', style='bar', colour=colors.green, altcolour=colors.red)

point_every=100
windowsize=10000
countposition=0
dataset=[]
#totstring=getTotalString(origin_gb_file_dir+orgabb,largest_contig,"CDS")
number_of_points=len(parsed_gb)/point_every
#while countposition<len(parsed_gb):
#    dataset.append((countposition,getWindowScore(totalstring,countposition,windowsize)))
#    countposition+=point_every





size=(len(parsed_gb.seq)/float(5400000))*26

gd_diagram.draw(format="circular", circular=True, pagesize=(size*cm,size*cm), start=0, end=200000, circle_core=0.5,tracklines=False)
#gd_diagram.draw(format="circular", circular=True, pagesize=(size*cm,size*cm), start=0, end=len(parsed_gb), circle_core=0.5,tracklines=False)
#gd_diagram.draw(format="circular", circular=True, pagesize=(20*cm,20*cm), start=0, end=len(parsed_gb))
gd_diagram.write("/Users/security/science/genomedraw.pdf", "PDF")
#gd_diagram.write("plasmid_linear.eps", "EPS")
#gd_diagram.write("plasmid_linear.svg", "SVG")
#gd_diagram.write("plasmid_linear.png", "PNG")




#SeqIO.write(set, output_handle, "genbank")

print "done"

