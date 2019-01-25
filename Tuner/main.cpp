#include "mbed.h"
#include "math.h"
#define pi 3.14

float samples[2048] = {0};

Serial pc(USBTX, USBRX);
//AnalogIn ulaz(PTB0);

DigitalOut Red(LED1);
DigitalOut Green(LED2);
DigitalOut Blue(LED3);
/*
int index(float* sum){

for(k = 0; k < 2048 / 8; k++)
{
    max_value2 = 0 ;
    //downsample 3 times by /1, /2, and /3
    float sum = magFFT[k] * magFFT[2*k] * magFFT[3*k];
    // find fundamental frequency by finding largest peak value
    if( sum > max_value2 && k > 3 )
    {
        //max_value2 = sum[k];
        fund_freq = k;
    }
}
}
*/
void vFFT(float data[], unsigned int nn)
{
    unsigned int n,mmax,m,j,istep,i;
    double wtemp,wr,wpr,wpi,wi,theta;
    float tempr,tempi;

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

    n=nn << 1;
    j=1;
    for (i=1; i<n; i+=2) {
        if(j>i) {
            SWAP(data[j],data[i]);
            SWAP(data[j+1],data[i+1]);
        }
        m=n >> 1;
        while (m >= 2 &&j>m) {
            j-=m;
            m >>= 1;
        }
        j+=m;
    }

    mmax=2;
    while (n > mmax) {
        istep=mmax << 1;
        theta=(6.28318530717959/mmax);
        wtemp=sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi=sin(theta);
        wr=1.0;
        wi=0.0;
        for (m=1; m<mmax; m+=2) {
            for (i=m; i<=n; i+=istep) {
                j=i+mmax;
                tempr=wr*data[j]-wi*data[j+1];
                tempi=wr*data[j+1]+wi*data[j];
                data[j]=data[i]-tempr;
                data[j+1]=data[i+1]-tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr=(wtemp=wr)*wpr-wi*wpi+wr;
            wi=wi*wpr+wtemp*wpi+wi;
        }
        mmax=istep;
    }

}

float getHammingValue(int i,int size)
{
    return (float) (0.54 - 0.46 * cos( (2*pi * i)/(size -1))) ;
}

void applyWindow(float* sample, int size)
{
    for(int i(0); i<size; i++) {
        sample[i] = sample[i] * getHammingValue(i,size);
    }
}

void Indicate(double frekf)
{
    const float tolerancija = 1;
    const float E(84.53) , A(110.77), D(147.8), G(195.31), B(244.86), e(325.03);
    const float EA(96.205) ,AD(128.4), DG(171.4), GB(221.45), Be(288.25);

    if(frekf <= E-tolerancija) {
        //low E
        pc.printf("E2-\n");
        Red = 0;
        Green = 1;
        Blue = 1;
    } else if(frekf > E-tolerancija && frekf <= E+tolerancija) {
        //E
        pc.printf("E2\n");
        Red = 1;
        Green = 0;
        Blue = 1;
    } else if(frekf>E+tolerancija && frekf <=EA ) {
        //high E
        pc.printf("E2+\n");
        Red = 1;
        Green = 1;
        Blue = 0;
    } else if(frekf > EA && frekf <= A-tolerancija) {
        //low A
        pc.printf("A2-\n");
        Red = 0;
        Green = 1;
        Blue = 1;
    } else if(frekf > A-tolerancija &&  frekf <= A+tolerancija) {
        //A
        pc.printf("A2\n");
        Red = 1;
        Green = 0;
        Blue = 1;
    } else if(frekf > A+tolerancija && frekf <= AD) {
        //high A
        pc.printf("A2+\n");
        Red = 1;
        Green = 1;
        Blue = 0;
    } else if(frekf > AD && frekf <= D-tolerancija) {
        //low D
        pc.printf("D3-\n");
        Red = 0;
        Green = 1;
        Blue = 1;
    } else if(frekf > D-tolerancija &&  frekf <= D+tolerancija) {
        //D
        pc.printf("D3\n");
        Red = 1;
        Green = 0;
        Blue = 1;
    } else if(frekf > D+tolerancija && frekf <= DG) {
        //high D
        pc.printf("D3+\n");
        Red = 1;
        Green = 1;
        Blue = 0;
    } else if(frekf > DG && frekf <= G-tolerancija) {
        //low G
        pc.printf("G3-\n");
        Red = 0;
        Green = 1;
        Blue = 1;
    } else if(frekf > G-tolerancija &&  frekf <= G+tolerancija) {
        //G
        pc.printf("G3\n");
        Red = 1;
        Green = 0;
        Blue = 1;
    } else if(frekf > G+tolerancija && frekf <= GB) {
        //high G
        pc.printf("G3+\n");
        Red = 1;
        Green = 1;
        Blue = 0;
    } else if(frekf > GB && frekf <= B-tolerancija) {
        //low B
        pc.printf("B3-\n");
        Red = 0;
        Green = 1;
        Blue = 1;
    } else if(frekf > B-tolerancija &&  frekf <= B+tolerancija) {
        //B
        pc.printf("B3\n");
        Red = 1;
        Green = 0;
        Blue = 1;
    } else if(frekf > B+tolerancija && frekf <= Be) {
        //high B
        pc.printf("B3+\n");
        Red = 1;
        Green = 1;
        Blue = 0;
    } else if(frekf > Be && frekf <= e-tolerancija) {
        //low e
        pc.printf("E4-\n");
        Red = 0;
        Green = 1;
        Blue = 1;
    } else if(frekf > e-tolerancija &&  frekf <= e+tolerancija) {
        //e
        pc.printf("E4\n");
        Red = 1;
        Green = 0;
        Blue = 1;
    } else {
        //high e
        pc.printf("E4+\n");
        Red = 1;
        Green = 1;
        Blue = 0;
    }
}

int main()
{
    pc.printf("Ready\n");
    wait(1);
    //while(1) {
        AnalogIn ulaz(PTB3);
        for(int i = 0; i<2048; i++) {
            samples[i] = 0;
        }

        for(int i = 0; i<2048; i++) {
            samples[i] = (float)ulaz;
            wait(0.000125);
        }
        //pc.printf("\nDone\n");

        applyWindow(samples, 2048);
        vFFT(samples, 1024);

        int index= 0 ;

        for(int i(0); i < 2048; i++)
            if(samples[i]< 0)samples[i]*=-1 ;

        for(int i(3) ; i < 2048; i++)
            if(samples[index]< samples[i])index = i ;

        float frekf = index * 8000./2048;
        //pc.printf("\n%f\n",frekf*0.373126);
        Indicate(frekf*0.373126);

        wait(1);
   // }
}