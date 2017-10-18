/**
  This header file is added to enable the barest minimum for organizing the code. This file is used 
  to modify the phase function. The user should return a variable theta.
**/

//Henyey-Greenstein Phase Function, "Handbook of Optical 
                       //Biomedical Diagnostics",2002,Chap3,p234, also see Boas2002

                       if(tmp0>EPS){  //if prop.g is too small, the distribution of theta is bad
		           tmp0=(1.f-prop.g*prop.g)/(1.f-prop.g+2.f*prop.g*rand_next_zangle(t));
		           tmp0*=tmp0;
		           tmp0=(1.f+prop.g*prop.g-tmp0)/(2.f*prop.g);

                           // when ran=1, CUDA gives me 1.000002 for tmp0 which produces nan later
                           // detected by Ocelot,thanks to Greg Diamos,see http://bit.ly/cR2NMP
                           tmp0=fmax(-1.f, fmin(1.f, tmp0));

		           theta=acosf(tmp0);
		           stheta=sinf(theta);
		           ctheta=tmp0;
                       }else{
			   theta=acosf(2.f*rand_next_zangle(t)-1.f);
                           sincosf(theta,&stheta,&ctheta);
                       }