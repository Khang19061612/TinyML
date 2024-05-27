#include <Arduino.h>
#include "Filter.h"
float omega0;
float dt;
int order =2;
float Q;
float domega;
float tn1 = 0;
float y_n[] = {0,0,0,0,0};
float x_n[] ={0,0,0,0};
float a_low[2];
float b_low[3];
float a_high[2];
float b_high[3];
float a_band[3];
float b_band[3];
float a_stop[3];
float b_stop[3];
//for moving average*********************************************
//***************************************************************
void lowpass_init(float f0, float fs)
{
  // f0: cutoff frequency (Hz)
  // fs: sample frequency (Hz)
  omega0 = 6.28318530718*f0;
  dt = 1.0/fs;
  tn1 = -dt;
  for(int k = 0; k < order+1; k++){
        x_n[k] = 0;
        y_n[k] = 0;        
      }
  float alpha = omega0*dt;
  float alphaSq = alpha*alpha;
  float beta[] = {1, sqrt(2), 1};
  float D = alphaSq*beta[0] + 2*alpha*beta[1] + 4*beta[2];
  b_low[0] = alphaSq/D;
  b_low[1] = 2*b_low[0];
  b_low[2] = b_low[0];
  a_low[0] = -(2*alphaSq*beta[0] - 8*beta[2])/D;
  a_low[1] = -(beta[0]*alphaSq - 2*beta[1]*alpha + 4*beta[2])/D;
}
void highpass_init(float f0, float fs)
{
  // f0: cutoff frequency (Hz)
  // fs: sample frequency (Hz)
  omega0 = 6.28318530718*f0;
  dt = 1.0/fs;
  tn1 = -dt;
  for(int k = 0; k < order+1; k++){
        x_n[k] = 0;
        y_n[k] = 0;        
      }
  float alpha = omega0*dt;
  float dtSq = dt*dt;
  float c[] = {omega0*omega0, sqrt(2)*omega0, 1};
  float D = c[0]*dtSq + 2*c[1]*dt + 4*c[2];
  b_high[0] = 4.0/D;
  b_high[1] = -8.0/D;
  b_high[2] = 4.0/D;
  a_high[0] = -(2*c[0]*dtSq - 8*c[2])/D;
  a_high[1] = -(c[0]*dtSq - 2*c[1]*dt + 4*c[2])/D;   
}
void bandpass_init(float f0, float fw, float fs)
{
  // f0: central frequency (Hz)
  // fw: bandpass width (Hz)
  // fs: sample frequency (Hz)
  omega0 = 6.28318530718*f0;
  domega = 6.28318530718*(fw+00);
  Q = omega0/domega;
  dt = 1.0/fs;
  tn1 = -dt;
  for(int k = 0; k < order+1; k++){
        x_n[k] = 0;
        y_n[k] = 0;    
        a_band[k]=0;
        b_band[k]=0;    
      }
  float alpha = omega0*dt;
  float D = pow(alpha,2) + 2*alpha/Q + 4;
  b_band[0] = 2*alpha/(Q*D);
  b_band[1] = 0;
  b_band[2] = -b_band[0];
  a_band[0] = 0;
  a_band[1] = -(2*pow(alpha,2) - 8)/D;
  a_band[2] = -(pow(alpha,2) - 2*alpha/Q + 4)/D;    
}
void bandstop_init(float f0, float fw, float fs)
{
  // f0: central frequency (Hz)
  // fw: bandpass width (Hz)
  // fs: sample frequency (Hz)
  omega0 = 6.28318530718*f0;
  domega = 6.28318530718*(fw+00);
  Q = omega0/domega;
  dt = 1.0/fs;
  tn1 = -dt;
  for(int k = 0; k < order+1; k++){
        x_n[k] = 0;
        y_n[k] = 0;    
        a_stop[k]=0;
        b_stop[k]=0;    
      }
  float alpha = omega0*dt;
  float D = pow(alpha,2) + 2*alpha/Q + 4;
  b_stop[0] = (pow(alpha,2)+4)/D;
  b_stop[1] = (2*pow(alpha,2)-8)/D;
  b_stop[2] = b_stop[0];
  a_stop[0] = 0;
  a_stop[1] = -(2*pow(alpha,2) - 8)/D;
  a_stop[2] = -(pow(alpha,2) - 2*alpha/Q + 4)/D;    
}
int lowpass(int x1)
{
  y_n[2]=y_n[1];
  y_n[1]=y_n[0];
  y_n[0]=0;
  y_n[0] = a_low[0]*y_n[1] +  a_low[1]*y_n[2] + b_low[0]*x1 + b_low[1]*x_n[0] + b_low[2]*x_n[1];
  x_n[2]=x_n[1];
  x_n[1]=x_n[0];
  x_n[0]=x1;
  return y_n[0];
}
int highpass(int x1)
{
  y_n[2]=y_n[1];
  y_n[1]=y_n[0];
  y_n[0]=0;
  y_n[0] = a_high[0]*y_n[1] +  a_high[1]*y_n[2] + b_high[0]*x1 + b_high[1]*x_n[0] + b_high[2]*x_n[1];
  x_n[2]=x_n[1];
  x_n[1]=x_n[0];
  x_n[0]=x1;
  return y_n[0];
}

int bandpass(int x1)
{

  y_n[2]=y_n[1];
  y_n[1]=y_n[0];
  y_n[0] = a_band[0]*y_n[0] +  a_band[1]*y_n[1] +  a_band[2]*y_n[2] +  b_band[0]*x1 + b_band[1]*x_n[0] + b_band[2]*x_n[1] ;
  x_n[1]=x_n[0];
  x_n[0]=x1;
  return y_n[0];
}

int bandstop(int x1)
{

  y_n[2]=y_n[1];
  y_n[1]=y_n[0];
  y_n[0] = a_stop[0]*y_n[0] +  a_stop[1]*y_n[1] +  a_stop[2]*y_n[2] +  b_stop[0]*x1 + b_stop[1]*x_n[0] + b_stop[2]*x_n[1] ;
  x_n[1]=x_n[0];
  x_n[0]=x1;
  return y_n[0];
}
int movingAvg(int *ptrArrNumbers, long *ptrSum, int pos, int len, int nextNum)
{
  //Subtract the oldest number from the prev sum, add the new number
  *ptrSum = *ptrSum - ptrArrNumbers[pos] + nextNum;
  //Assign the nextNum to the position in the array
  ptrArrNumbers[pos] = nextNum;
  //return the average
  return *ptrSum / len;
}

void reset_filter()
{
	y_n[0]=y_n[1]=y_n[2]=y_n[3]=y_n[4]=0;
	x_n[0]=x_n[1]=x_n[2]=x_n[3]=0;
}
