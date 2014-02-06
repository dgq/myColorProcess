#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "time.h"

#define N 3
#define D8 255
#define mGamma 2.2
#define CHUNK 1024

static const float mColorXform[3][3] = 
{
	{ 0.2990f,  0.5870f,  0.1140f},
	{-0.1690f, -0.3310f,  0.5000f},
	{ 0.5000f, -0.4190f, -0.0810f}
};

static float mCCMatrix[3][3] =
{
	{ 1.7537f, -0.5628f, -0.1909f},
	{-0.2188f,  1.3846f, -0.1658f},
	{ 0.0214f, -0.5475f,  1.5260f}
};

static float mCCMatrixInvert[3][3]=
{
	{0.604870f,	0.288163f,  0.106977f},	
	{0.098814f,	0.801727f,  0.099469f},	
	{0.026970f,	0.283604f,  0.689495f}
};

static float mRGBBeforeGamma[3][3]=
{
	{186.00f,  70.00f, 73.00f},  	//Red block
	{88.00f, 160.00f, 87.00f},  	//Green block
	{41.33f, 15.56f, 16.22f}  		//Red block

};


static float mRGBBeforeCCM[3][3]=
{
	{186.00f,  70.00f, 73.00f},  	//Red block
	{88.00f, 160.00f, 87.00f},  	//Green block
	{41.33f, 15.56f, 16.22f}  		//Red block

};

static float mRGBAfterGamma[3][3]=
{
	{186.00f,  70.00f, 73.00f},  	//Red block
	{88.00f, 160.00f, 87.00f},  	//Green block
	{65.00f, 75.00f, 163.00f}  		//Red block
};

static float mHSVA[3][3]=    		//Hue:0-360 S:0-1   V:0-1  RGB After gamma
{
	{186.00f,  70.00f, 73.00f},  	//Red block
	{88.00f, 160.00f, 87.00f},  	//Green block
	{41.33f, 15.56f, 16.22f}  		//Red block
};

static float mHSVB[3][3]=    		//Hue:0-360 S:0-1   V:0-1 RGB Before gamma
{
	{186.00f,  70.00f, 73.00f},  	//Red block
	{88.00f, 160.00f, 87.00f},  	//Green block
	{41.33f, 15.56f, 16.22f}  		//Red block
};


void mFlush(void)
{
	int c;
	while ( (c = getchar()) != '\n' && c != EOF ) ;
}

/*************************Matrix invert begin*******************************/

float MatDet(float *p, int n)
{
    int r, c, m;
    int lop = 0;
    float result = 0;
    float mid = 1;

    if (n != 1)
    {
        lop = (n == 2) ? 1 : n;           
        for (m = 0; m < lop; m++)
        {
            mid = 1;            
            for (r = 0, c = m; r < n; r++, c++)
            {
                mid = mid * (*(p+r*n+c%n));
            }
            result += mid;
        }
        for (m = 0; m < lop; m++)
        {
            mid = 1;           
            for (r = 0, c = n-1-m+n; r < n; r++, c--)
            {
                mid = mid * (*(p+r*n+c%n));
            }
            result -= mid;
        }
    }
    else
        result = *p;
    return result;
}

float Creat_M(float *p, int m, int n, int k)
{
    int len;
    int i, j;
    float mid_result = 0;
    int sign = 1;
    float *p_creat, *p_mid;

    len = (k-1)*(k-1);            
    p_creat = (float*)malloc(len*sizeof(float)); 
    p_mid = p_creat;
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < k; j++)
        {
            if (i != m && j != n) 
            {
                *p_mid++ = *(p+i*k+j);
            }
        }
    }
    sign = (m+n)%2 == 0 ? 1 : -1;    
    mid_result = (float)sign*MatDet(p_creat, k-1);
    free(p_creat);
    return mid_result;
}

void print(float *p, int n)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        
        for (j = 0; j < n; j++)
        {
            printf("%.3f\t",*p++);
        }
        printf("\n");
    }
}

int mGauss(float A[][N], float B[][N], int n)
{
    int i, j, k;
    float mmax, temp;
    float t[N][N];                
    
    for (i = 0; i < n; i++)        
    {
        for (j = 0; j < n; j++)
        {
            t[i][j] = A[i][j];
        }
    }
    
    for (i = 0; i < n; i++)        
    {
        for (j = 0; j < n; j++)
        {
            B[i][j] = (i == j) ? (float)1 : 0;
        }
    }
    for (i = 0; i < n; i++)
    {
        
        mmax = t[i][i];
        k = i;
        for (j = i+1; j < n; j++)
        {
            if (fabs(t[j][i]) > fabs(mmax))
            {
                mmax = t[j][i];
                k = j;
            }
        }
        
        if (k != i)
        {
            for (j = 0; j < n; j++)
            {
                temp = t[i][j];
                t[i][j] = t[k][j];
                t[k][j] = temp;
                
                temp = B[i][j];
                B[i][j] = B[k][j];
                B[k][j] = temp;
            }
        }
        
        if (t[i][i] == 0)
        {
            printf("\nThere is no inverse matrix!\n");
            return 1;
        }
        
        temp = t[i][i];
        for (j = 0; j < n; j++)
        {
            t[i][j] = t[i][j] / temp;        
            B[i][j] = B[i][j] / temp;        
        }
        for (j = 0; j < n; j++)        
        {
            if (j != i)                
            {
                temp = t[j][i];
                for (k = 0; k < n; k++)        
                {
                    t[j][k] = t[j][k] - t[i][k]*temp;
                    B[j][k] = B[j][k] - B[i][k]*temp;
                }
            }
        }
    }
    return 0;
}

int mVerifyMatrix(float aa[][N], float bb[][N], int n)
{
	int i, j, k;
	float a[n][n], b[n][n], c[n][n];
	float mSum = 0;

	for (i=0; i<n; i++)
	{
		for (j=0; j<n; j++)
		{
			a[i][j] = aa[i][j];
			b[i][j] = bb[i][j];
		}
	}
	#if 0
	//for (i=0; i<n; i++)
	{
		
			c[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0];
			c[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1] + a[0][2]*b[2][1];
			c[0][2] = a[0][0]*b[0][2] + a[0][1]*b[1][2] + a[0][2]*b[2][2];

			c[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0] + a[1][2]*b[2][0];
			c[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1] + a[1][2]*b[2][1];
			c[1][2] = a[1][0]*b[0][2] + a[1][1]*b[1][2] + a[1][2]*b[2][2];

			c[2][0] = a[2][0]*b[0][0] + a[2][1]*b[1][0] + a[2][2]*b[2][0];
			c[2][1] = a[2][0]*b[0][1] + a[2][1]*b[1][1] + a[2][2]*b[2][1];
			c[2][2] = a[2][0]*b[0][2] + a[2][1]*b[1][2] + a[2][2]*b[2][2];

	}
	#else
	
	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
		{
			mSum=0;
			for(k=0;k<n;k++)
				mSum+=a[i][k]*b[k][j];
			c[i][j]=mSum;
		}

    #endif

	printf("Verify matix invert:\n" );
	for (i=0; i<n; i++)
	{
		for (j=0; j<n; j++)
			printf("%10.6f\t", c[i][j] );
		printf("\n");
	}
	for (i=0; i<n; i++)
	{
		if (fabs(c[i][i] - 1.000000f)<=1e-6)
		{
			printf("Matrix invert operation is Correct, success!!!\n");
			return 0;
			
		}
		else
		{
			printf("Matrix invert operation is incorrect!\n");
			return 1;
		}
	}
}

int mArrayInvert(void)
{
	#if 0
	float *buffer, *p;            
    int row, num;                
    int i, j;
    float determ;                
    float a[N][N], b[N][N];
    int n;
	int mRet = 0;
    
	printf("Test invert matrix using invert matrix defination!\n");
	printf("Please input matrix Row number:\t");
	//mFlush();
	scanf("%d", &row);
	//mFlush();
    num = 2 * row * row;
    buffer = (float *)malloc(num*sizeof(float));        
    p = buffer;
    if (NULL != p)
    {
        for (i = 0; i < row; i++)
        {
			printf("Please input the value of %d row: \t", i+1);
            for (j = 0; j < row; j++)
            {
				//mFlush();
				scanf("%f", p++);
				//mFlush();
            }
        }
    }
    else
    {
	   printf("Can't distribute memory\n");
    }
    printf("The original matrix : \n");
    print(buffer, row);                
    
    determ = MatDet(buffer, row);    
    p = buffer + row * row;
    if (determ != 0)
    {
		printf("The determinant of the matrix is ", determ);
        for (i = 0; i < row; i++)    
        {
            for (j = 0; j < row; j++)
            {
                *(p+j*row+i) = Creat_M(buffer, i, j, row)/determ;
            }
        }
        printf("The inverse matrix is:\n"); 
        print(p, row);                
    }
    else
    {
        printf("The determinant is 0, and there is no inverse matrix!\n");
    }
    free(buffer);   
	
	#endif

	          
                   
    int i, j;
    float determ;                
    float a[N][N], b[N][N];
    int n = 3;
	int mRet = 0;     
	printf("\nTo get convert matrix using mGauss!\n");
   
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
			//scanf("%f", &a[i][j]);
			a[i][j] = mCCMatrix[i][j];
        }
    }
	printf("Orignal Matrix:\n");
	for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
			printf("%10.6f\t", a[i][j]);
        }
		printf("\n");
    }
   
    if (!mGauss(a, b, n))
    {
		printf("Inverse matrix:\n");
        for (i = 0; i < n; i++)
        {
            
            for (j = 0; j < n; j++)
            {
                printf("%10.6f\t", b[i][j]);
            }
           printf("\n");
        }

		printf("Copy Inverse Matrix:\n");
		for (i=0; i<n; i++)
		{
			for (j=0; j<n; j++)
			{
				mCCMatrixInvert[i][j] = b[i][j];
				printf("%10.6f\t", mCCMatrixInvert[i][j]);
			}
			printf("\n");
		}
    } 

	mRet = mVerifyMatrix(a, b, n);
	
	return mRet;
}

/*************************Matrix invert end*********************************/

/***********************gamma correction begin******************************/
void mGammaCorrection(void)
{
	int i, j, m;
    float p[3] = {0.0f, 0.0f, 0.0f};
	printf("\nInput Decoding Gamma(r > 1):  %3.2f\n", mGamma);
	//scanf("%f", &mGamma);
	printf("Output Encoding Gamma(r< 1):  %3.2f\n\n", 1.0f/mGamma);

	
	for (m=0; m<N; m++)
	{
		switch (m)
		{
		case 0:
			printf("Red Block:\n");
			break;
		case 1:
			printf("Green Block:\n");
			break;
		case 2:
			printf("Blue Block:\n");
			break;
		default:
			printf("Unknown Block\n");
			break;
		}
		for (i = 0; i < N; i++) 
		{
			//scanf("%f", p[i]);
			p[i] = mRGBBeforeGamma[m][i];
			printf("%10.3ff,", p[i]);
		}
		printf("\t(0-255(before gamma)\n");
	
		for (i = 0; i < N; i++) 
		{
			p[i]=p[i]/D8;
			if (p[i]/D8<=0.018f)
				p[i]=4.5*p[i];
			else
				p[i]=1.099*(pow(p[i],1.0f/mGamma));
			if (p[i]>=1.0f) 
				p[i]=1.0f;
			p[i] = p[i] * D8; 
			mRGBAfterGamma[m][i] = p[i];
			printf("%10.3ff,", mRGBAfterGamma[m][i]);
		}
		printf("\t(0-255)(after gamma)\n");
		
	}
}

void mDeGamma(void)
{
	int i, m;
	float p[3] = {0.0f, 0.0f, 0.0f};
	printf("\n************************************************\n");
	printf("\nInput Decoding Gamma(r > 1):  %3.2f\n", mGamma);
	//scanf("%f", &mGamma);
	printf("Output Encoding Gamma(r< 1):  %3.2f\n\n", 1.0f/mGamma);

	for (m=0; m<N; m++)
	{
		switch (m)
		{
		case 0:
			printf("Red Block:\n");
			break;
		case 1:
			printf("Green Block:\n");
			break;
		case 2:
			printf("Blue Block:\n");
			break;
		default:
			printf("Unknown Block\n");
			break;
		}
		for (i = 0; i < N; i++) 
		{
			//scanf("%f", p[i]);
			p[i] = mRGBAfterGamma[m][i];
			printf("%10.3ff,", p[i]);
		}
		printf("\t(0-255)(after gamma)\n");

		
		for (i = 0; i < N; i++) 
		{
			p[i]=p[i]/D8;
			if (p[i]/D8 <= 4.5f*0.018f)
				p[i]= p[i]/4.5f;
			else
				p[i]=exp(log(p[i]/1.099)*mGamma);  //  1/mGamma
			if (p[i]>=1.0f) 
				p[i]=1.0f;
			p[i] = p[i] * D8; 
			mRGBBeforeGamma[m][i] = p[i];
			printf("%10.3ff,", mRGBBeforeGamma[m][i]);
		}
		printf("\t(0-255)(before gamma)\n");
	}
	printf("\n************************************************\n");
	
}

/*************************gamma correction end******************************/

/****************************RGB->HSV begin********************************/
void mRGB2HSV(void)
{
	int i, m;
	float mR, mG, mB;
	float mMax, mMin, mDelta;
	float mHue, mSat, mValue;
	float mTemp;

	for (m=0; m<2*N; m++)
	{
		switch (m)
		{
		case 0:
			printf("Red block(before gamma):\n");
			printf("RGB Value:");
			mR = mRGBBeforeGamma[0][0]; 
			mG = mRGBBeforeGamma[0][1];
			mB = mRGBBeforeGamma[0][2];
			printf("%10.3ff,%10.3ff,%10.3ff\n", mR, mG, mB);
			break;
		case 1:
			printf("Red block(after gamma):\n");
			printf("RGB Value:");
			mR = mRGBAfterGamma[0][0]; 
			mG = mRGBAfterGamma[0][1];
			mB = mRGBAfterGamma[0][2];
			printf("%10.3ff,%10.3ff,%10.3ff\n", mR, mG, mB);
			break;
		case 2:
			printf("Green block(before gamma):\n");
			printf("RGB Value:");
			mR = mRGBBeforeGamma[1][0]; 
			mG = mRGBBeforeGamma[1][1];
			mB = mRGBBeforeGamma[1][2];
			printf("%10.3ff,%10.3ff,%10.3ff\n", mR, mG, mB);
			break;
		case 3:
			printf("Green block(after gamma):\n");
			printf("RGB Value:");
			mR = mRGBAfterGamma[1][0]; 
			mG = mRGBAfterGamma[1][1];
			mB = mRGBAfterGamma[1][2];
			printf("%10.3ff,%10.3ff,%10.3ff\n", mR, mG, mB);
			break;
		case 4:
			printf("Blue block(before gamma):\n");
			printf("RGB Value:");
			mR = mRGBBeforeGamma[2][0]; 
			mG = mRGBBeforeGamma[2][1];
			mB = mRGBBeforeGamma[2][2];
			printf("%10.3ff,%10.3ff,%10.3ff\n", mR, mG, mB);
			break; 
		case 5:
			printf("Blue block(after gamma):\n");
			printf("RGB Value:");
			mR = mRGBAfterGamma[2][0]; 
			mG = mRGBAfterGamma[2][1];
			mB = mRGBAfterGamma[2][2];
			printf("%10.3ff,%10.3ff,%10.3ff\n", mR, mG, mB);
			break; 
		default:
			printf("Unknown mRGB2HSV RGB Data!\n");
			break;
		}
		
		mR = mR/D8; mG = mG/D8; mB = mB/D8;

		if (mR>=mG)
		{
			mMax = mR;
			mMin = mG;
		}
		else
		{
			mMax = mG;
			mMin = mR;
		}
		if (mMax <= mB) 
			mMax = mB;
		if (mMin >= mB)
			mMin = mB;

	//Hue calculation

		if (mMax == mMin)
			mHue = 0.0f;
		else
		{
			mTemp = 60.0f/(mMax - mMin);
			if ((mMax == mR)&&(mG>=mB))
				mHue = mTemp*(mG-mB)+0.0f;
			if ((mMax==mR)&&(mG<mB))
				mHue = mTemp*(mG-mB)+360.0f;
			if (mMax==mG)
				mHue = mTemp*(mB-mR)+120.0f;
			if (mMax==mB)
				mHue = mTemp*(mR-mG)+240.0f;
		}
	//Saturation Calculation

		if (mMax==0.0f)
			mSat = 0.0f;
		else
			mSat = 1 - mMin/mMax;

	//Value Calculation

		mValue = mMax;
		
	//Print HSV value
		switch (m)
		{
		case 0:
			printf("HSV Value:");
			mHSVB[0][0]= mHue; 
			mHSVB[0][1]= mSat; 
			mHSVB[0][2]= mValue;
			printf("%10.3ff,%10.3ff,%10.3ff\n", mHue, mSat, mValue);
			break;
		case 1:
			printf("HSV Value:");
			mHSVA[0][0]= mHue; 
			mHSVA[0][1]= mSat; 
			mHSVA[0][2]= mValue;
			printf("%10.3ff,%10.3ff,%10.3ff\n\n", mHue, mSat, mValue);
			break;
		case 2:
			printf("HSV Value:");
			mHSVB[1][0]= mHue; 
			mHSVB[1][1]= mSat; 
			mHSVB[1][2]= mValue;
			printf("%10.3ff,%10.3ff,%10.3ff\n", mHue, mSat, mValue);
			break;
		case 3:
			printf("HSV Value:");
			mHSVA[1][0]= mHue; 
			mHSVA[1][1]= mSat; 
			mHSVA[1][2]= mValue;
			printf("%10.3ff,%10.3ff,%10.3ff\n\n", mHue, mSat, mValue);
			break;
		case 4:
			printf("HSV Value:");
			mHSVB[2][0]= mHue; 
			mHSVB[2][1]= mSat; 
			mHSVB[2][2]= mValue;
			printf("%10.3ff,%10.3ff,%10.3ff\n", mHue, mSat, mValue);
			break;
		case 5:
			printf("HSV Value:");
			mHSVA[2][0]= mHue; 
			mHSVA[2][1]= mSat; 
			mHSVA[2][2]= mValue;
			printf("%10.3ff,%10.3ff,%10.3ff\n\n", mHue, mSat, mValue);
			break; 
		default:
			printf("Unknown mRGB2HSV HSV Data!\n");
			break;
		}
	}
}

/****************************RGB->HSV begin********************************/

/****************************HSV->RGB begin********************************/
void mHSV2RGB()
{
	int i, m;
	int mH; 
	float mF, mP, mQ, mT;
	float mHue, mSat, mValue;
	float mR, mG, mB;

	
	//Note: H:[0-360), S:[0-1), V:[0-1)
	for (m=0; m<2*N; m++)
	{
		switch (m)
		{
		case 0:
			printf("\nRed block(before gamma):\n");
			printf("HSV Value:");
			mHue = mHSVB[0][0]; 
			mSat = mHSVB[0][1]; 
			mValue = mHSVB[0][2];
			printf("%10.3ff,%10.3ff,%10.3ff\n", mHue, mSat, mValue);
			break;
		case 1:
			printf("Red block(after gamma):\n");
			printf("HSV Value:");
			mHue = mHSVA[0][0]; 
			mSat = mHSVA[0][1]; 
			mValue = mHSVA[0][2];
			printf("%10.3ff,%10.3ff,%10.3ff\n", mHue, mSat, mValue);
			break;
		case 2:
			printf("\nGreen block(before gamma):\n");
			printf("HSV Value:");
			mHue = mHSVB[1][0]; 
			mSat = mHSVB[1][1]; 
			mValue = mHSVB[1][2];
			printf("%10.3ff,%10.3ff,%10.3ff\n", mHue, mSat, mValue);
			break;
		case 3:
			printf("Green block(after gamma):\n");
			printf("HSV Value:");
			mHue = mHSVA[1][0]; 
			mSat = mHSVA[1][1]; 
			mValue = mHSVA[1][2];
			printf("%10.3ff,%10.3ff,%10.3ff\n", mHue, mSat, mValue);
			break;
		case 4:
			printf("\nBlue block(before gamma):\n");
			printf("HSV Value:");
			mHue = mHSVB[2][0]; 
			mSat = mHSVB[2][1]; 
			mValue = mHSVB[2][2];
			printf("%10.3ff,%10.3ff,%10.3ff\n", mHue, mSat, mValue);
			break;
		case 5:
			printf("Blue block(after gamma):\n");
			printf("HSV Value:");
			mHue = mHSVA[2][0]; 
			mSat = mHSVA[2][1]; 
			mValue = mHSVA[2][2];
			printf("%10.3ff,%10.3ff,%10.3ff\n", mHue, mSat, mValue);
			break; 
		default:
			printf("Unknown mRGB2HSV HSV Data!\n");
			break;
		}

		mH=(int)mHue/60;
		mF=mHue/60.0f-(float)mH;
		mP=mValue*(1-mSat);
		mQ=mValue*(1-mF*mSat);
		mT=mValue*(1-(1-mF)*mSat);

		if (mH==0)
		{
			mR=mValue; mG=mT; mB=mP;
		}
		if (mH==1)
		{
			mR=mQ; mG=mValue; mB=mP;
		}
		if (mH==2)
		{
			mR=mP; mG=mValue; mB=mT;
		}
		if (mH==3)
		{
			mR=mP; mG=mQ; mB=mValue;
		}
		if (mH==4)
		{
			mR=mT; mG=mP; mB=mValue;
		}
		if (mH==5)
		{
			mR=mValue; mG=mP; mB=mQ;
		}

		mR = mR*D8; mG = mG*D8; mB = mB*D8;
		switch (m)
		{
		case 0:
			printf("RGB Value:");
			mRGBBeforeGamma[0][0] = mR; 
			mRGBBeforeGamma[0][1] = mG;
			mRGBBeforeGamma[0][2] = mB;
			printf("%10.3ff,%10.3ff,%10.3ff\n", mR, mG, mB);
			break;
		case 1:
			printf("RGB Value:");
			mRGBAfterGamma[0][0] = mR; 
			mRGBAfterGamma[0][1] = mG;
			mRGBAfterGamma[0][2] = mB;
			printf("%10.3ff,%10.3ff,%10.3ff\n", mR, mG, mB);
			break;
		case 2:
			printf("RGB Value:");
			mRGBBeforeGamma[1][0] = mR; 
			mRGBBeforeGamma[1][1] = mG;
			mRGBBeforeGamma[1][2] = mB;
			printf("%10.3ff,%10.3ff,%10.3ff\n", mR, mG, mB);
			break;
		case 3:
			printf("RGB Value:");
			mRGBAfterGamma[1][0] = mR; 
			mRGBAfterGamma[1][1] = mG;
			mRGBAfterGamma[1][2] = mB;
			printf("%10.3ff,%10.3ff,%10.3ff\n", mR, mG, mB);
			break;
		case 4:
			printf("RGB Value:");
			mRGBBeforeGamma[2][0] = mR; 
			mRGBBeforeGamma[2][1] = mG;
			mRGBBeforeGamma[2][2] = mB;
			printf("%10.3ff,%10.3ff,%10.3ff\n", mR, mG, mB);
			break; 
		case 5:
			printf("RGB Value:");
			mRGBAfterGamma[2][0] = mR; 
			mRGBAfterGamma[2][1] = mG;
			mRGBAfterGamma[2][2] = mB;
			printf("%10.3ff,%10.3ff,%10.3ff\n", mR, mG, mB);
			break; 
		default:
			printf("Unknown mRGB2HSV RGB Data!\n");
			break;
		}
	}
}

/*****************************HSV->RGB end*********************************/

/***********************Input data from File Begin**************************/
int mInputDataFromFile(void)
{
	FILE* pF;
	int mSize = 0;
	int mFileSize = 0;
	int mRet = 0;
	char* mBuffer=NULL;
	char mFileName[] = "myColorData.txt";
	

	char mValidChar[13] = 
	{'-','f','.','0','1','2','3','4','5','6','7','8','9'};
	char mEndofData = '~';
	char mTemp;

	pF = fopen(mFileName, "rb");

	if (pF==NULL) 
	{
		mRet = 1;
		goto mError; 
	}
	else
	{
		printf("File %s Open Success!\n", mFileName);
		
		printf("\n************************************************\n");
		do
		{
			mTemp = fgetc(pF);
			printf("%c", mTemp);
			if (mTemp == mEndofData)
				mSize = mFileSize + 1;
			if (mTemp != EOF)
				mFileSize++; 
			
			if (mSize > CHUNK)
			{
				mRet = 4;
				goto mError;
			}
		}while(mTemp!=EOF);
		printf("\n************************************************\n");
		printf("\nFile size %d, actual size %d\n", mFileSize, mSize);
		if (mSize == 0)
		{
			mRet = 2;
			goto mError; 
		}
		else
			mBuffer = (char *)malloc(sizeof(char)*mSize); 
		rewind(pF);
		mRet=fread(mBuffer, 1, mSize, pF);
		if (mRet!=mSize)
		{
			printf("\nfread mRet %d mSize %d\n", mRet, mSize);
			mRet = 3;
			goto mError; 
		}


	}
	
	free(mBuffer);
	mBuffer = NULL;
	printf("Buffer %s, %p\n", mBuffer, mBuffer);
	
	fclose(pF);
	printf("File Point %s, %p\n", pF, pF);
	printf("\n************************************************\n");
	return 0;
	mError:
	{
		if (mBuffer != NULL) 
		{
			printf("\nSomething wrong for Free Buffer, try again!\n");
			free(mBuffer);
			printf("Buffer %s, %p\n", mBuffer, mBuffer);
		}
			
		if (pF != NULL)
		{
			printf("\nSomething wrong for close File, try again!\n");
			fclose(pF);
			printf("File Point %s, %p\n", pF, pF);
		}
		switch (mRet)
		{
		case 1:
			printf("Error open file!\n");
			break;
		case 2:
			printf("File size incorrect!\n");
			break;
		case 3:
			printf("File size mismatch!\n");
			break;
		case 4:
			printf("File Get Char error!\n");
			break;
		default:
			printf("unknown error!\n");
			break;
		}
		return 1;
		
	}
	
}
/***********************Input data from File End****************************/

/***********************Get RGB data Before CCM E***************************/
void mGetRGBBeforeCCM(void)
{
	int i, m;
	float p[3] = {0.0f, 0.0f, 0.0f};

	printf("\n************************************************\n");
	
	for (m=0; m<N; m++)
	{
		switch (m) 
			{
			case 0:
				printf("\nRed Block:\n");
				break;
			case 1:
				printf("\nGreen Block:\n");
				break;
			case 2:
				printf("\nBlue Block:\n");
				break;
			default:
				printf("Unknown Block\n");
				break;
			}
			for (i = 0; i < N; i++) 
			{
				p[i] = mRGBBeforeGamma[m][i]/D8;
				printf("%10.3ff,", mRGBBeforeGamma[m][i]);
			}
			printf("\t(0-255)(Before gamma)\n");

			for (i = 0; i < N; i++) 
			{
				mRGBBeforeCCM[m][i] =  mCCMatrixInvert[i][0]*p[0]
									  +mCCMatrixInvert[i][1]*p[1]
									  +mCCMatrixInvert[i][2]*p[2];
				
				mRGBBeforeCCM[m][i] = mRGBBeforeCCM[m][i] * D8;
				printf("%10.3ff,", mRGBBeforeCCM[m][i]);
					 
			}
			printf("\t(0-255)(Before CCM)\n");
	}
	printf("\n************************************************\n");
	
}

/***********************Get RGB data Before CCM X***************************/


/*****************************Print to File E*******************************/
int mPrint2File(void)
{
	int i,j;
	char mFileName[30];
	FILE* mFp;
	time_t mRAWTime;
	struct tm* mCurrentTime;

	time(&mRAWTime);
	mCurrentTime =  localtime(&mRAWTime);
	strftime(mFileName, sizeof(mFileName), "my_%Y%m%d%H%M%S.txt",mCurrentTime);
	//It's strange, if use "%D", it will show nothing!

	mFp=fopen(mFileName, "wb");
	if (mFp==NULL)
	{
		printf("Open file error!");
		return 1;
	}
	fprintf(mFp,"%s\n",asctime(mCurrentTime));
	fprintf(mFp, "Only for test\n");

	fclose(mFp);

	return 0;
}


/*****************************Print to File X*******************************/

int main(void)
{
	int mRet = 0;
	int mOp = 0;
	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	printf("~                                                            ~\n");
	printf("~                                                            ~\n");
	printf("~  Color Proecss Function                                    ~\n");
	printf("~  Editor: dgq                                               ~\n");
	printf("~  Version 0.1        Date: Feb 3, 2014                      ~\n");
	printf("~  Hello, listed below is code to process color, enjoy it!!! ~\n");
	printf("~                                                            ~\n");
	printf("~                                                            ~\n");
	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

	//Read the color data from the text file named myColorData.txt
	#if 0
	mRet = mInputDataFromFile();
	if (mRet!=0)
	{
		printf("mInputDataFromFile error!\n");
		return 0;
	}
	#endif

	do
	{
		printf("\n\n\n0: Test;\n");
		printf("1: Calculate inverse of CCMatrix;\n");
		printf("2: Calculate RGB value before gamma from the value after;\n");
		printf("3: Calculate RGB value after gamma from the value before;\n");
		printf("4: Calculate HSV value from RGB value;\n");
		printf("5: Calculate RGB value from HSV value;\n");
		printf("6: Calculate RGB value Before CCM;\n\n");
		printf("\nInput function:  "); 
		scanf("%d", &mOp);
		
		
		switch (mOp) 
		{
		case 0:
			printf("**********************Test E**************************\n");
			mPrint2File();
			printf("**********************Test X**************************\n");
			break;
		case 1:
			printf("*****************Inverse Matrix E*********************\n");
			mRet = mArrayInvert();
			if (mRet != 0)
			{
				printf("mArray Invert error!\n");
				return 1; 
			}
			printf("*****************Inverse Matrix X*********************\n");
			break;
		case 2:
			printf("*************RGB(gamma before<-after) E***************\n");

			mDeGamma();

			printf("*************RGB(gamma before<-after) X***************\n");
			break;
		case 3:
			printf("*************RGB(gamma after<-before) E***************\n");

			mGammaCorrection();

			printf("*************RGB(gamma after<-before) X***************\n");
			break;
		case 4:
			printf("********************RGB->HSV E************************\n");

			mRGB2HSV();

			printf("********************RGB->HSV X************************\n");
			break;
		case 5:
			printf("********************HSV->RGB E************************\n");

			mHSV2RGB();

			printf("********************HSV->RGB X************************\n");
			break;
		case 6:
			printf("****************RGB value before CCM******************\n");

			mGetRGBBeforeCCM();

			printf("****************RGB value before CCM******************\n");
			break;
		default:
			printf("***************To be continued loading****************\n");
			break;
		}

	} while ((mOp>=0)&&(mOp<=6));
	return 0;
}

