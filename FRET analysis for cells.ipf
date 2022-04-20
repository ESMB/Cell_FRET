#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
macro run()

multiple()

end

function multiple()
variable numberoffolders=(3)
make/o/n=(numberoffolders) pixels,fretpixels,centre,width,fraction
make/o/n=(numberoffolders)/T filelist
string filename="red.tif"      
string filename2="acceptor.tif"      
string filename3="donor.tif"

string root="Macintosh HD:Users:Mathew:Dropbox (Cambridge University):FRET for Mathew:190920 Figure1_Eiii:6days:A30P:"

variable a
for(a=0;a<(numberoffolders);a+=1)
filelist[a]=root+num2str(a)+":"
endfor

variable f
for(f=0;f<numberoffolders;f+=1)
kill()
	setdatafolder root:
	string folder=num2str(f)
	//PRINT FOLDER
	
	//Load all of the images
	newdatafolder/s $folder
	string path=filelist[f]
	string toload=path+filename
	
	print toload
	
	ImageLoad/T=tiff/Q/N=red  toload
	string toload2=path+filename2
	ImageLoad/T=tiff/Q/N=acceptor  toload2
	
	string toload3=path+filename3
	ImageLoad/T=tiff/Q/N=donor  toload3
	
	analyse_images(path)

wave width,centre,number_of_fret_pixels,number_of_red_pixels,fraction
variable wi=width[0]
variable ce=centre[0]
variable nufre=number_of_fret_pixels[0]
variable nure=number_of_red_pixels[0]
variable fr=fraction[0]

wave binary_FRET
newimage binary_FRET
ModifyImage binary_FRET ctab= {*,*,Grays,0}
ModifyGraph tick=3
ModifyGraph noLabel=2
string save1=filelist[f]+"binary_FRET"
SavePICT/O/E=-2 as save1


wave FRET_image
newimage FRET_image
ModifyImage FRET_image ctab= {*,*,YellowHot,0}
ModifyImage FRET_image ctab= {0,1,YellowHot,0}
ModifyGraph tick=3
ModifyGraph noLabel=2
ColorScale/C/N=text0/F=0/A=RC/E image=FRET_image
string savefret=filelist[f]+"FRET_Image"
SavePICT/O/E=-2 as savefret

wave Z_image
newimage Z_image
ModifyImage Z_image ctab= {*,*,Rainbow,0}
ModifyGraph tick=3
ModifyGraph noLabel=2
ModifyImage Z_image ctab= {-2,2,RedWhiteBlue,0}
ColorScale/C/N=text0/F=0/A=RC/E image=Z_image
string savez=filelist[f]+"Z_Image"
SavePICT/O/E=-2 as savez

wave binary_red
newimage binary_red
ModifyImage binary_red ctab= {*,*,Grays,0}
ModifyGraph tick=3
ModifyGraph noLabel=2
string save2=filelist[f]+"binary_red"
SavePICT/O/E=-2 as save2



setdatafolder root:
wave pixels,fretpixels,centre,width,fraction
pixels[f]=nure
fretpixels[f]=nufre
centre[f]=ce
width[f]=wi
fraction[f]=fr

endfor
string tosave3=root+"data.txt"
Save/J/M="\n"/W/DLIM=","/O filelist,fraction,fretpixels,pixels,width,centre as tosave3
end




function analyse_images(path)
string path
wave red,acceptor,donor
// Thresholds to extract
variable red_threshold=200		// Threshold red channel, since this should be independent of FRET
variable acceptor_threshold=0	// This is the threshold that should be applied to decide whether there is a signal in the FRET acceptor channel
variable a,b,c
make/o/n=(512,512) thresholded_red=0,binary_red=0,binary_FRET=0,FRET_image=0,Z_image=0,S_image=0
make/o/n=1 xcoord,ycoord,redI,donorI,acceptorI,FRET,zparam,Stoich
make/o/n=1 fretcount=0
for(a=0;a<(512);a+=1)
	for(b=0;b<512;b+=1)
		if(red[a][b]>red_threshold)
			redimension/n=(c+1) xcoord,ycoord,redI,donorI,acceptorI,FRET,zparam,stoich
			xcoord[c]=a
			ycoord[c]=b
			redI[c]=red[a][b]
			donorI[c]=donor[a][b]
			acceptorI[c]=acceptor[a][b]
			FRET[c]=acceptor[a][b]/(donor[a][b]+acceptor[a][b])
			zparam[c]=ln(acceptor[a][b]/donor[a][b])
			stoich[c]=(donor[a][b]+acceptor[a][b])/(donor[a][b]+acceptor[a][b]+red[a][b])
			if(acceptorI[c]>acceptor_threshold)
			fretcount+=1
			binary_fret[a][b]=1
			FRET_image[a][b]=FRET[c]
			Z_image[a][b]=zparam[c]
			S_image[a][b]=stoich[c]
			endif
			c+=1
			thresholded_red[a][b]=red[a][b]
			binary_red[a][b]=1
			
		endif
	endfor
endfor
Make/N=40/O zparam_Hist;DelayUpdate
Histogram/B={-4,0.2,40} zparam,zparam_Hist;DelayUpdate
Display zparam_Hist
ModifyGraph mode=5
Label left "# Events";DelayUpdate
Label bottom "Z"
CurveFit gauss zparam_Hist /D 
string saveas=path+"Z_Hist"
SavePICT/O/E=-2 as saveas
wave w_coef
make/o/n=1 centre=w_coef[2],width=w_coef[3]
make/o/n=1 number_of_red_pixels=c,number_of_fret_pixels=fretcount
make/o/n=1 fraction=fretcount/c

Make/N=20/O FRET_Hist;DelayUpdate
Histogram/B={0,0.05,20} FRET,FRET_Hist;DelayUpdate
Display FRET_Hist
Label left "# Events";DelayUpdate
Label bottom "FRET Efficiency"
ModifyGraph mode=5
string saveas2=path+"FRET_Hist"
SavePICT/O/E=-2 as saveas2

string FRETpic=path+"FRET_32bit.tif"

ImageSave/T="tiff"/DS=32 FRET_image as FRETpic

end



/////////////////////////////

/////////////////////////////

/////////////////////////////

/////////////////////////////
function analyse_afterwards()
variable eccenthresh=0.5
variable lengththresh=50
variable A_param=0.000201077
variable b_param=0.59178009
variable c_param=11.16478573
variable d_param=0.220696
setdatafolder root:
wave filelist
variable l=dimsize(filelist,0)

make/o/n=(l) oligomers,fibrils,oligomers_conc,fibrils_conc,fibrils_by_mass,oligomers_by_mass

variable a,b,c

for(a=0;a<(l);a+=1)
	string folder=num2str(a)
	setdatafolder $folder
	make/o/n=(1) fibril_lengths,oligomer_lengths
	wave eccen,length
	variable oligomercount=0,fibrilcount=0
	
	for(b=0;b<(dimsize(eccen,0));b+=1)
		if(eccen[b]<eccenthresh | length[b]>lengththresh)
			fibrilcount+=1
			redimension/n=(fibrilcount+1) fibril_lengths
			fibril_lengths[fibrilcount]=length[b]
		else
		
			oligomercount+=1
			redimension/n=(oligomercount+1) oligomer_lengths
			oligomer_lengths[oligomercount]=length[b]
		endif
	endfor
	Make/N=200/O length_Hist;DelayUpdate
Histogram/B={0,20,200} length,length_Hist;DelayUpdate

Make/N=200/O oligomer_length_Hist;DelayUpdate
Histogram/B={0,20,200} oligomer_lengths,oligomer_length_Hist;DelayUpdate

Make/N=200/O fibril_length_Hist;DelayUpdate
Histogram/B={0,20,200} fibril_lengths,fibril_length_Hist;DelayUpdate

make/o/n=200 length_hist_norm
wavestats length_hist
variable h
for(h=0;h<200;h+=1)
length_hist_norm[h]=length_hist[h]/v_sum
endfor

make/o/n=1 test
make/o/n=200 oligomer_length_hist_norm
wavestats oligomer_length_hist
for(h=0;h<200;h+=1)
oligomer_length_hist_norm[h]=oligomer_length_hist[h]/v_sum
endfor

make/o/n=200 fibril_length_hist_norm
wavestats fibril_length_hist
for(h=0;h<200;h+=1)
fibril_length_hist_norm[h]=fibril_length_hist[h]/v_sum
endfor

// Need to do the mass part now- take bin-centre value as the length. So 10 nm, 30 nm, 50 nm...
// Need to multiply length by the number of monomers per nm - one monomer = 0.47 nm fibril length- so length/0.47 	

variable g
make/o/n=200 oligomer_norm_mass
for(g=0;g<(200);g+=1)
oligomer_norm_mass[g]=oligomer_length_hist_norm[g]*((g)*20+10)/0.47		// This is number of monomers
endfor

make/o/n=200 fibril_norm_mass
for(g=0;g<(200);g+=1)
fibril_norm_mass[g]=fibril_length_hist_norm[g]*((g)*20+10)/0.47				// This is number of monomers
endfor

make/o/n=1 fibril_avg,oligomer_avg

wavestats fibril_norm_mass

fibril_avg[0]=v_sum					// Average monomer units (weighted)
variable fibrilav=fibril_avg[0]
wavestats oligomer_norm_mass		// // Average monomer units (weighted)
oligomer_avg[0]=v_sum
variable oligav=oligomer_avg[0]

setdatafolder root:

oligomers[a]=oligomercount			
fibrils[a]=fibrilcount

oligomers_conc[a]=C_param*((A_param-D_param)/((oligomers[a]/(67.3^2))-D_param)-1)^(1/B_param)

fibrils_conc[a]=C_param*((A_param-D_param)/((fibrils[a]/(67.3^2))-D_param)-1)^(1/B_param)

fibrils_by_mass[a]=fibrils_conc[a]*fibrilav
oligomers_by_mass[a]=oligomers_conc[a]*oligav



endfor

			
			
	

end


function anal2()
setdatafolder root:
wave filelist
variable l=dimsize(filelist,0)

make/o/n=(l) lengthsskel,lengthloc,monomerunits,monomerunitsfromlength

variable a,b,c

for(a=0;a<(l);a+=1)
	string folder=num2str(a)
	setdatafolder $folder
	howmanymonomers()
	
	
	wave monomers,lengthfromlocs,length
	
	wavestats monomers
	
	variable mon=v_avg
	 wavestats lengthfromlocs
	 variable lenlo=v_avg
	 wavestats length
	 variable len=v_avg
	 make/o/n=(dimsize(length,0)) monomerlen
	for(b=0;b<(dimsize(length,0));b+=1)
		monomerlen[b]=length[b]/0.47
	endfor
	
	wavestats monomerlen
	variable monlen=v_avg
	setdatafolder root:
	
	lengthsskel[a]=len
	lengthloc[a]=lenlo
	monomerunits[a]=mon
	monomerunitsfromlength[a]=monlen
endfor	

end


function howmanymonomers()
variable numberofframes=4000
// These are from the linear fit
variable aparam=0
variable bparam=1.4456 

// Name of the wave containing localisations:

wave pointspercluster

duplicate/o pointspercluster,locs


// Gubbins here:

variable a,b,c

make/o/n=(dimsize(locs,0)) lengthfromlocs,monomers

for(a=0;a<(dimsize(locs,0));a+=1)

lengthfromlocs[a]=aparam+bparam*(4000/numberofframes)*locs[a]		// From fit
monomers[a]=	lengthfromlocs[a]/0.47														// 1 monomer = 0.47 nm of length. 

endfor



end


function anal()
wave wave0,wave1

duplicate/o wave0,length
duplicate/o wave1,eccen

make/o/n=(20,40) matrix=0

variable a,b,c

variable lengthsep=100
variable eccensep=0.025
variable lengthlow=0,ecclow=0
make/o/n=20 lengthaxis
make/o/n=40 eccaxis
for(a=0;a<20;a+=1)
	
	variable lengthhigh=lengthlow+lengthsep
	for(b=0;b<40;b+=1)
		variable ecchigh=ecclow+eccensep
			for(c=0;c<(dimsize(length,0));c+=1)
				if(eccen[c]>ecclow && eccen[c]<ecchigh && length[c]<lengthhigh && length[c]>lengthlow)
				matrix[a][b]+=1
				endif
			endfor
			eccaxis[b]=(ecclow+ecchigh)/2
			ecclow+=eccensep
			
	//print ecclow
	endfor
	lengthaxis[a]=(lengthlow+lengthhigh)/2
lengthlow+=lengthsep
ecclow=0
//print lengthlow
endfor
	Display;DelayUpdate
AppendMatrixContour matrix vs {lengthaxis,eccaxis}

ModifyContour matrix labels=0,autoLevels={*,*,100}
ModifyContour matrix ctabLines={*,*,Rainbow,1}
ModifyGraph width=283.465,height=283.465
ModifyGraph mirror=1,minor=1,notation=1
Label bottom "Length (nm)"
Label left "Eccentricity"
ModifyGraph gmSize=20
SetAxis bottom 50,2000
SetAxis left 0,1
ModifyGraph gfSize=20,gmSize=0
ModifyGraph btLen=5,stLen=3



end






function kill()

	variable                      winMask;
 
	variable                      i,n;
	variable                      all=0x1000+0x40+0x10+0x4+0x2+0x1;
	string                        theWins;
 
	winMask = !winMask ? all : winMask;
 
	theWins = winList("*",";","WIN:"+num2iStr(winMask & all));
	for(i=0,n=itemsInList(theWins,";") ; i<n ; i+=1)
		doWindow/K $stringFromList(i,theWins,";");
	endfor;
end

macro close_windows()
kill()
end


function ext()
wave filelist
make/o/n=(200,(dimsize(filelist,0))) fiblength,oliglength
setdatafolder root:

variable l=dimsize(filelist,0)


variable a,b,c

for(a=0;a<(l);a+=1)
	string folder=num2str(a)
	setdatafolder $folder
wave fibril_length_hist,oligomer_length_hist

duplicate/o fibril_length_hist,:fibtemp
duplicate/o oligomer_length_hist,:oligomertemp
setdatafolder root:
wave oligomertemp,fibtemp
for(c=0;c<200;c+=1)
fiblength[c][a]=fibtemp[c]
oliglength[c][a]=oligomertemp[c]
endfor


endfor

		setdatafolder root:	
			
	

end






function all_lengths_verses_locs()
setdatafolder root:
wave filelist
variable l=dimsize(filelist,0)

make/o/n=1 clusterlength,clustercount

variable a,b,c,le

for(a=0;a<(l);a+=1)
	string folder=num2str(a)
	setdatafolder $folder
wave eccen
wavestats/q eccen
// First of all, need to determine number of localisations per cluster
print (dimsize(eccen,0))
	wave clusternumber,length

	wavestats/q clusternumber
	
	variable numberclust=v_max
	
	make/o/n=(numberclust) pointspercluster
	
	
	Histogram/B={1,1,numberclust} Clusternumber,pointspercluster
	
	for(c=0;c<(dimsize(pointspercluster,0));c+=1)
		variable pc=pointspercluster[c]
		variable lc=length[c]
	setdatafolder root:	
		redimension/n=(le+1) clusterlength,clustercount
		clusterlength[le]=lc
		clustercount[le]=pc
		le+=1
	setdatafolder $folder
	endfor
	
	setdatafolder root:


endfor

			
			
	

end





function moreaccuratelength()
wave clustercount,clusterlength

variable thresholdlength=30
variable thresholdcounts=100
make/o/n=1 thresh_counts,thresh_length

variable a,b,c

for(a=0;a<(dimsize(clustercount,0));a+=1)
	if(clusterlength[a]>thresholdlength&&clustercount[a]>thresholdcounts)
		redimension/n=(b+1) thresh_counts,thresh_length
		
		thresh_counts[b]=clustercount[a]
		thresh_length[b]=clusterlength[a]
		
		b+=1
	endif
endfor


end



function allmeans()
wave oligomers,fibrils,oligomers_conc,fibrils_conc,fibrils_by_mass,oligomers_by_mass,lengths,monomerunitsfromlength

duplicate/o oligomers,input
means()
wave replicamean,replicastdev
duplicate/o replicamean,oligomers_mean
duplicate/o replicastdev,oligomers_stdev

duplicate/o fibrils,input
means()
wave replicamean,replicastdev
duplicate/o replicamean,fibrils_mean
duplicate/o replicastdev,fibrils_stdev

duplicate/o oligomers_conc,input
means()
wave replicamean,replicastdev
duplicate/o replicamean,oligomers_conc_mean
duplicate/o replicastdev,oligomers_conc_stdev

duplicate/o fibrils_conc,input
means()
wave replicamean,replicastdev
duplicate/o replicamean,fibrils_conc_mean
duplicate/o replicastdev,fibrils_conc_stdev

duplicate/o oligomers_by_mass,input
means()
wave replicamean,replicastdev
duplicate/o replicamean,oligomers_by_mass_mean
duplicate/o replicastdev,oligomers_by_mass_stdev

duplicate/o fibrils_by_mass,input
means()
wave replicamean,replicastdev
duplicate/o replicamean,fibrils_by_mass_mean
duplicate/o replicastdev,fibrils_by_mass_stdev

duplicate/o lengths,input
means()
wave replicamean,replicastdev
duplicate/o replicamean,lengths_mean
duplicate/o replicastdev,lengths_stdev

duplicate/o monomerunitsfromlength,input
means()
wave replicamean,replicastdev
duplicate/o replicamean,monomerunitsfromlength_mean
duplicate/o replicastdev,monomerunitsfromlength_stdev





end 
function means()
variable repeats=4
variable replica=3

// Input wave

wave input

//duplicate/o fibrils_by_mass,input

variable a,b,c

make/o/n=(dimsize(input,0)/4) meansfromrepeat

for(a=0;a<(dimsize(input,0));a+=repeats)
	variable add=0
	for(b=0;b<repeats;b+=1)
		add=add+input[a+b]
		print input[a+b]
	endfor
	print 0
	meansfromrepeat[c]=add/repeats
c+=1
endfor
c=0
make/o/n=((dimsize(meansfromrepeat,0)/replica)) replicamean,replicastdev

for(a=0;a<(dimsize(meansfromrepeat,0));a+=replica)
make/o/n=(replica) temp
	for(b=0;b<replica;b+=1)
		temp[b]=meansfromrepeat[a+b]
	endfor
wavestats temp
string name=num2str(c)
duplicate/o temp,$name

replicamean[c]=v_avg
replicastdev[c]=v_sdev
c+=1
endfor

end


function timeaxismake()
make/o/n=8 timeaxis={0,12,23,3,6,9} 



end



Function twogauss(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = y0 + a1*exp(-((x-xc1)/w1)^2) + a2*exp(-((x-xc2)/w2)^2)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 7
	//CurveFitDialog/ w[0] = y0
	//CurveFitDialog/ w[1] = a1
	//CurveFitDialog/ w[2] = xc1
	//CurveFitDialog/ w[3] = w1
	//CurveFitDialog/ w[4] = a2
	//CurveFitDialog/ w[5] = xc2
	//CurveFitDialog/ w[6] = w2

	return w[0] + w[1]*exp(-((x-w[2])/w[3])^2) + w[4]*exp(-((x-w[5])/w[6])^2)
End
