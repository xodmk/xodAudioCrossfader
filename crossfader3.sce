// /////////////////////////////////////////////////////////////////////////////////////////////
// #############################################################################################
// begin : header
// #############################################################################################
// /////////////////////////////////////////////////////////////////////////////////////////////
//
// __<<name=> "odmkXfade1.sce">>__

// ___::((JIROBATA Programming Industries))::___
// ___::((ODMK:odorousbeast:BarutanBreaks:djoto:2014:2015:2016))::___
// ___::((created by eschei))___

// impllementation of a STFT - short-time fourier transform <-> iSTFT re-synthesis
// rotate bins between each frame


// include components:

// ___::((equal power crossfade))::___
// ___::((24 but .wav file read & write))::___
//

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// This script implements crossfading between two .wav files
// reads in two .wav files, crossfades, writes out mixed/faded .wav file
// allows for contour control
// uses look-up table with linear interpolation to generate attenuation curves
// Multiple plots are created to show contour curves for multiple values of n

// Real time controls:
// n_pos - array of contour values - 1 = cos/-cos, 2...n increases curve sharpness
// x_pos - crossfade position - corresponds to sample position in .wav file
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// /////////////////////////////////////////////////////////////////////////////////////////////
// #############################################################################################
// end : header
// #############################################################################################
// /////////////////////////////////////////////////////////////////////////////////////////////

xdel(winsid()) //-> closes all open graphs (= matlab.close('all'))
clear
clc


exec('/Users/apple/odmk-djoto/odmk-sci/odmk_code/scilab/wav_rdwr_scilab/odmkWav_rw.sce');


///////////////////////////////////////////////////////////////////////////////////////////////
//begin : stack adjustment
///////////////////////////////////////////////////////////////////////////////////////////////

//try for ~10MBytes
//stacksize input is doubles
//
max_bytes = 10*10e6;
max_bits=16;
//bits=24;
bytespersample=ceil(max_bits/8);
max_data_bytes=max_bytes-(12+24);    //total size - header,format,etc. approx
max_stack=max_data_bytes/8;

//if size > max then
//    error('wav file too large');
//else
stacksize(max_stack)
//    stacksize('max')
//    stacksize('min')
//    sz=stacksize()
//end

///////////////////////////////////////////////////////////////////////////////////////////////
//end : stack adjustment
///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//begin : function definitions
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//begin : table functions
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

////////////////////////////////////////////////////////////////
//function: tablegen
//generate look-up table

function [tbout] = tablegen(depth)
    
    //create look-up table entries for different waveforms
    //define an n deep table
    //depth=n;
    tl=(0:1/(depth-1):1); 
    // store 1 cycle 
    for q=1:length(tl)
        table1(q)=cos(2*%pi/4*(tl(q)));    //2 multiple necessary because function range changes from [0:2] to [0:1] for normalization
    end
    tbout = table1;
    

endfunction


////////////////////////////////////////////////////////////////
//function: xtable
//look-up table using linear interpolation
//assumes a table 'tb' has been generated and initialized with correct function

function [iout] = xtable(addr,td)
    //assume the input addr is a value betwee 0+ and 1  
    //assumes that the table depth is < than 2^dataWidth  
    //table depth = td
    //qntWidth = ceil(log2(td));    
    
    
    //assume address input is a value between 0 and 1 - scale to log2(depth) address bits 
    
    
    xaddr=addr*(2^(log2(td))-2)+1;   //scale addr to match table depth -> +1 for scilab 1 - n addressing 
    qaddr=floor(xaddr);                //quantize for use as table address
    //look up table outputs
    //tb is generated from a call from Main to tablegen so the table isn't rebuilt each time
    yHigh=tb(qaddr+1);
    yLow=tb(qaddr);
    
    //linear interpolation
    //p(x) = f(x0) + (f(x1)-f(x0))/(x1-x0)*(x-x0)    //x1-x0 = 1
    //p(x) = f(x0) + (f(x1)-f(x0))*(x-x0)
    //iout = y0 + (y1-y0)*(xaddr-saddr);
    
    //***check this, dependent on difference between table depth and address sweep range?????
    
    if (qaddr==td) then    //avoid addr overflow for linear interpolation addressing - should investigate
        iout = tb(qaddr);
    else
        yHigh=tb(qaddr+1);
        yLow=tb(qaddr);
        //linear interpolation
        ////p(x) = f(x0) + (f(x1)-f(x0))/(x1-x0)*(x-x0)    //x1-x0 = 1
        ////p(x) = f(x0) + (f(x1)-f(x0))*(x-x0)
        iout = yLow + (yHigh-yLow)*(xaddr-qaddr);
    end

endfunction
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//end : table functions
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


////////////////////////////////////////////////////////////////
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//begin : crossfade functions
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

////////////////////////////////////////////////////////////////
//function: xfade1
//crossfader function - direct implementation
function [fx, fy, cpower] = xfade1(x, n)
    
    fx=0;    //n rows by x columns
    fy=0;
    cpower=0;

    fx=cos(%pi/4*((2*x-1)^(2*n-1)+1));
    fy=cos(%pi/4*((1-2*x)^(2*n-1)+1));
    //constant power equation:
    cpower=fx^2+fy^2;

endfunction
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////
//function: xfade2
//crossfader function - interpolated look-up table implementation
//inputs: x = position (-1<->1), n = contour value (range?), td = LUT depth  
function [fx2, fy2, cpower2] = xfade2(x, n, td)

    
    //initialize internals
    intx1=0;
    inty1=0;
    intn1=0;
    intx2=0;
    inty2=0;
    
    //initialize outputs
    fx2=0;
    fy2=0;
    cpower2=0;
    
    intx1=2*x-1;
    inty1=1-2*x;
    intn1=2*n-1;
    
    //multiply-accumulator used to implement exponential
    //consider serial vs. parallel implementations for hardware implementation
    
    //when x is normalized to [0:1], macc output (intx2,inty2) has range [0:2]
    for i=1:intn1
        if intn1==1 then
            intx2=intx1+1;
            inty2=inty1+1;
        elseif (intn1>1 & i==1) then
            intx2=intx1;
            inty2=inty1;
        elseif (intn1>1 & i==intn1) then
            intx2=intx2*intx1+1;
            inty2=inty2*inty1+1;
        else
            intx2=intx2*intx1;
            inty2=inty2*inty1;
        end
    end

    //temp, replace with look-up table
    //fx2_tap1=cos(%pi/4*(intx2));
    //fy2_tap1=cos(%pi/4*(inty2));
    
   
    //scale input
    //sintx2=intx2/2;
    //sinty2=inty2/2;
    
    //fx2=xtable(intx2/2,table_depth);    //scale inputs from range [0:2] -> [0:1]
    //fy2=xtable(inty2/2,table_depth);
    fx2=xtable(intx2/2,td);    //scale inputs from range [0:2] -> [0:1]
    fy2=xtable(inty2/2,td);    
    

    //constant power equation:
    cpower2=fx2^2+fy2^2;
    

endfunction
////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//end : crossfade functions
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//end : function definitions
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//begin : main script
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//begin : read .wav file
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//read existing audio file
//siz=xwaveread(the-lorax-tkno01.wav,'size');
wavfile_inA='/Users/apple/odmk-djoto/odmk-sci/odmk_code/scilab/phasevocoder/wavSource/601_90ALTOSAX_C_chop.wav';
wavfile_inB='/Users/apple/odmk-djoto/odmk-sci/odmk_code/scilab/phasevocoder/wavSource/The_Amen_Break_odmk.wav';


//[y,fs,bits]=xwavread("/odmk-djoto/odmk-dsp/scilab/wav_rdwr_scilab/atatakai_16b");
[wavA,fs,bits]=xwavread(wavfile_inA);
[wavB,fs,bits]=xwavread(wavfile_inB);

szA=length(wavA(:,1));
szB=length(wavB(:,1));
nchanA=length(wavA(1,:));
nchanB=length(wavB(1,:));
mprintf(' <xwavread out> # of audio samples = %i x # of channels = %i\n',szA,nchanA);
mprintf(' <xwavread fs output> Sample Rate = %f (Hz)\n',fs);
mprintf(' <xwavread bits output>  bits per sample = %us\n\n',bits)

wavADataL=wavA(:,1);
wavADataR=wavA(:,2);
wavBDataL=wavB(:,1);
wavBDataR=wavB(:,2);

if szA <= szB then
    Dlength = length(wavA(:,1));
else
    Dlength = length(wavB(:,1));
end
chADataL = wavA((1:Dlength),1);
chADataR = wavA((1:Dlength),2);
chBDataL = wavB((1:Dlength),1);
chBDataR = wavB((1:Dlength),2);

wavAClone=[chADataL,chADataR];
wavBClone=[chBDataL,chBDataR];

////////////////////////////////////////////////////////////////
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//end : read .wav file
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


//define array's for contour (n_pos), and x axis range (x_pos)
//n=[1,2,3,4,5,6,7,8,9,10,11,12,13];
n_pos=(1:1:13)
//x_pos=(0:1/511:1);
x_pos=(0:1/1023:1);

//dummy signals
w1=rand(1,2^13);    //2^13=8192 2^14=16384 
w2=rand(1,2^13);    //length=16384 



//cosine function - linear equal power crossfade function
cosinex=zeros(1,length(x_pos));
for m=1:length(x_pos)
    cosinex(m)=cos(x_pos(m));
end

//g(x) - equation used for adjusting contour
govx=zeros(length(n_pos),length(x_pos));
for p=1:length(n_pos)
    for l=1:length(x_pos)
      govx(p,l)=%pi/4*(2*x_pos(l)-1)^(2*n_pos(p)-1)+1;
    end
end



y1=zeros(length(n_pos),length(x_pos));
y2=zeros(length(n_pos),length(x_pos));
cp1=zeros(length(n_pos),length(x_pos));

z1=zeros(length(n_pos),length(x_pos));
z2=zeros(length(n_pos),length(x_pos));

cp2=zeros(length(n_pos),length(x_pos));


//build table

//set Look-Up Table Depth for function look up
table_depth=4096;    //2^13
tb = tablegen(table_depth)

for i=1:length(n_pos)
    for j=1:length(x_pos)
	[y1(i,j),y2(i,j),cp1(i,j)] = xfade1(x_pos(j),n_pos(i));
     //[z1(i,j),z2(i,j),cp2(i,j)] = xfade2(x_pos(j),n_pos(i));
     [z1(i,j),z2(i,j),cp2(i,j)] = xfade2(x_pos(j),n_pos(i),table_depth);
    end
end

//check error between direct computation and look-up table 
error_int1=y1-z1;
error_int2=y2-z2;

//attenuate input signals w1,w2 & mix
//sweep fader position across length of signal
ts = (0:1/(length(chADataL)-1):1);
n_static=1;

//initialize crossfade arrays
xfA=zeros(1,length(chADataL));
xfB=zeros(1,length(chADataL));

chAL=zeros(1,length(chADataL));
chAR=zeros(1,length(chADataL));
chBL=zeros(1,length(chBDataL));
chBR=zeros(1,length(chBDataL));
mixA=zeros(1,length(chADataL));
mixB=zeros(1,length(chBDataL));

cptest=zeros(1,length(chADataL));

//
for s=1:length(chADataL)
    [xfA(s),xfB(s),cptest(s)] = xfade2(ts(s),n_static,table_depth);
    chAL(s) = chADataL(s)*xfA(s);
    chAR(s) = chADataR(s)*xfA(s);
    chBL(s) = chBDataL(s)*xfB(s);
    chBR(s) = chBDataR(s)*xfB(s);
    mixA(s) = chAL(s) + chBL(s);
    mixB(s) = chAR(s) + chBR(s);
end    

//calculate RMS
rms_chAL = sqrt(mean(chAL.^2));
rms_chAR = sqrt(mean(chAR.^2));
rms_chBL = sqrt(mean(chBL.^2));
rms_chBR = sqrt(mean(chBR.^2));
rms_mixA = sqrt(mean(mixA.^2));
rms_mixB = sqrt(mean(mixB.^2));


//format for writing .wav -> [n,2]
//mix outputs:
odmkMixWav = [mixA',mixB'];


///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//end : main script
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//begin : plotting
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////


//plotting
for k=1:length(n_pos)
    scf(1)
        plot(x_pos,y1(k,:),'red',x_pos,y2(k,:),'green',x_pos,z1(k,:),'blue',x_pos,z2(k,:),'magenta')
        title 'overlay of crossfader attenuation plots (ideal&LUT)) ; contour value n=[1-13]' 
        xlabel 'fader position (0=all ch1 <-mix-> 1=all ch2)', ylabel 'attenuation (blue=f(x) / magenta=f(1-x))'
    scf(2)
        plot(x_pos,z1(k,:),'cyan',x_pos,z2(k,:),'magenta')
        title 'interpolated table crossfader attenuation plots ; contour value n=[1-13]'
        xlabel 'fader position (0=all ch1 <-mix-> 1=all ch2)', ylabel 'attenuation (cyan=f(x) / majenta=f(1-x_pos))'
    scf(4)
        plot(x_pos,y1(k,:)-z1(k,:),'black')
        title 'error: direct computation vs interpolated look-up table ; contour value n=[1-13]'
        xlabel 'fader position (0=all ch1 <-mix-> 1=all ch2)', ylabel 'error (direct calculation vs interpolated table)'
//    scf(5)
//        plot(x_pos,error_int1(k,:),'black',x_pos,error_tb1(k,:),'red')
//        title 'error direct computation vs look-up table ; contour value n=[1-13]'
//        xlabel 'x position', ylabel 'error (black=y1-z1 / red=y1-z2_tap)'    
    scf(6)
        plot(x_pos,cp1(k,:),x_pos,cp2(k,:))
        title 'power plots: f(x)^2+f(1-x)^2   (overlay of 13 iterations)'
        xlabel 'fader position (0=all ch1 <-mix-> 1=all ch2)', ylabel 'power'
//    scf(7)
//        plot(x_pos,cosinex,'blue',x_pos,govx(k,:),'magenta')
//        title 'cos(x) and g(x) = pi/4*(2*x-1)^(2*n-1)+1;'
//        xlabel 'x position', ylabel 'blue=cos(x) / magenta=g(x)'

end


scf(8)
        plot(ts,chAL,'red',ts,chAR,'magenta')
        title 'attenuated signals:   ch1 output = red   ch2 output = magenta'
        xlabel 'fader position (0=all ch1 <-mix-> 1=all ch2)', ylabel 'attenuation'
scf(9)
        plot(ts,mixA,'red')
        title 'mix output signal   (overlay of 13 iterations)'
        xlabel 'fader position (0=all ch1 <-mix-> 1=all ch2) / mix output', ylabel 'attenuation'        
scf(10)
        plot(ts,cptest)
        title 'mix output power   (overlay of 13 iterations)'
        xlabel 'fader position (0=all ch1 <-mix-> 1=all ch2)', ylabel 'power (zoomed to 1)'  
        
        
///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//end : plotting
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//begin : write into wav files:
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////

//the-lorax-tkno01.wav

bits = 24;

//create a new audio file
odmkMixWav_out1='/Users/apple/odmk-djoto/odmk-sci/odmk_code/scilab/crossfader/wavaudio/odmk_xfmix1';
xwavwrite(odmkMixWav,fs,bits,odmkMixWav_out1);
mprintf(' <xwavwrite out> wrote .wav file %s\n\n',odmkMixWav_out1)

////create a new audio file
////nm='unicorn_clone1';    32bit IEEE float
//nm_out2='Users/apple/odmk-djoto/odmk-dsp/odmk_code/scilab/osc/odmk_sqrpulse';
//xwavwrite(???,fs,bits,nm_out2);
//mprintf(' <xwavwrite out> wrote .wav file %s\n\n',nm_out2)


///////////////////////////////////////////////////////////////////////////////////////////////
//#############################################################################################
//end : write into wav files:
//#############################################################################################
///////////////////////////////////////////////////////////////////////////////////////////////                 