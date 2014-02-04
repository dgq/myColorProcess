#include "stdio.h"
#include "math.h"
#include "stdlib.h"
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

static float mRGBAfterGamma[3][3]=
{
	{186.00f,  70.00f, 73.00f},  	//Red block
	{88.00f, 160.00f, 87.00f},  	//Green block
	{41.33f, 15.56f, 16.22f}  		//Red block
};

static float mHSV[3][3]=    		//Hue:0-360 S:0-1   V:0-1
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

bool Gauss(float A[][N], float B[][N], int n)
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
	printf("\nTo get convert matrix using Gauss!\n");
   
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
   
    if (!Gauss(a, b, n))
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
	int i, j;
    float p[3] = {0.0f, 0.0f, 0.0f};
	printf("\nInput Decoding Gamma(r > 1):  %3.2f\n", mGamma);
	//scanf("%f", &mGamma);
	printf("Output Encoding Gamma(r< 1):  %3.2f\n", 1.0f/mGamma);

	printf("\nRGB value before Gamma:\n");
	for (i=0; i<N; i++)
	{
		//scanf("%f", p[i]);
		p[i] = mRGBBeforeGamma[0][i];
		printf("%10.3f", p[i]);
	}
	printf("\t(0-255)\n");
	for (i=0; i<N; i++)
	{
		//scanf("%f", p[i]);
		mRGBBeforeGamma[1][i] =p[i]/D8;
		printf("%10.3f", mRGBBeforeGamma[1][i]);
	}
	printf("\t(0-1)\n");
	printf("\nRGB value after Gamma:\n");
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
		printf("%10.3f", p[i]);
	}
	printf("\t(0-255)\n");
	for (i=0; i<N; i++)
	{
		//scanf("%f", p[i]);
		mRGBAfterGamma[1][i] =p[i]/D8;
		printf("%10.3f", mRGBAfterGamma[1][i]);
	}
	printf("\t(0-1)\n");
}

void mDeGamma(void)
{
	int i;
	float p[3] = {0.0f, 0.0f, 0.0f};
	
	printf("\nInput Decoding Gamma(r > 1):  %3.2f\n", mGamma);
	//scanf("%f", &mGamma);
	printf("Output Encoding Gamma(r< 1):  %3.2f\n", 1.0f/mGamma);

	printf("\nRGB value After Gamma:\n");
	for (i=0; i<N; i++)
	{
		//scanf("%f", p[i]);
		p[i] = mRGBAfterGamma[0][i];
		printf("%10.3f", p[i]);
	}
	printf("\t(0-255)\n");
	for (i=0; i<N; i++)
	{
		//scanf("%f", p[i]);
		mRGBAfterGamma[1][i] =p[i]/D8;
		printf("%10.3f", mRGBAfterGamma[1][i]);
	}
	printf("\t(0-1)");
	

	printf("\nRGB value before Gamma:\n");
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
		printf("%10.3f", p[i]);
	}
	printf("\t(0-255)\n");
	for (i=0; i<N; i++)
	{
		//scanf("%f", p[i]);
		mRGBBeforeGamma[1][i] =p[i]/D8;
		printf("%10.3f", mRGBBeforeGamma[1][i]);
	}
	printf("\t(0-1)\n");
	
}

/*************************gamma correction end******************************/

/****************************RGB<->HSV begin********************************/
void mRGB2HSV(void)
{
	int i;
	float mR, mG, mB;
	float mMax, mMin, mDelta;
	float mHue, mSat, mValue;
	float mTemp;

	mR = mRGBBeforeGamma[1][0];
	mG = mRGBBeforeGamma[1][1];
	mB = mRGBBeforeGamma[1][2];

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

	mHSV[0][0]= mHue; mHSV[0][1]=mSat*100; mHSV[0][2]=mValue*100;
	mHSV[1][0]= mHue; mHSV[1][1]=mSat; mHSV[1][2]=mValue;
	
	printf("\nRGB value After Gamma:\n");
	for (i=0; i<N; i++)
		printf("%10.3f", mRGBBeforeGamma[0][i]);
	printf("\t(0-255)\n");
	for (i=0; i<N; i++)
		printf("%10.3f", mRGBBeforeGamma[1][i]);
	printf("\t(0-1)\n");
	printf("\nHSV value(H:0-360, S:0-100, V:0-100)\n");
	for (i=0; i<N; i++)
		printf("%10.3f", mHSV[0][i]);
	printf("\nHSV value(H:0-360, S:0-1, V:0-1)\n");
	for (i=0; i<N; i++)
		printf("%10.3f", mHSV[1][i]);
	printf("\n");
}

/****************************RGB<->HSV begin********************************/

/****************************HSV<->RGB begin********************************/
void mHSV2RGB()
{
	int i=0;
	int mH; 
	float mF, mP, mQ, mT;
	float mHue, mSat, mValue;
	float mR, mG, mB;

	//u can input the value directly.
	//Note: H:[0-360), S:[0-1), V:[0-1)
	mHue=mHSV[1][0]; mSat=mHSV[1][1]; mValue=mHSV[1][2];

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

	
	printf("\nHSV value(H:0-360, S:0-1, V:0-1)\n");
	for (i=0; i<N; i++)
		printf("%10.3f", mHSV[1][i]);

	printf("\n\nRGB value\n");
	printf("%10.3f", mR*255); printf("%10.3f", mG*255); 
	printf("%10.3f\t(0-255)\n", mB*255);
	printf("%10.3f", mR); printf("%10.3f", mG);
	printf("%10.3f\t(0-1)\n", mB);
	
}

/*****************************HSV<->RGB end*********************************/

/***********************Input data from File Begin**************************/
int mInputDataFromFile(void)
{
	FILE* pF;
	int mSize = 0;
	int mFileSize = 0;
	int mRet = 0;
	char* mBuffer=NULL;
	

	char mValidChar[13] = 
	{'-','f','.','0','1','2','3','4','5','6','7','8','9'};
	char mEndofData = '~';
	char mTemp;

	pF = fopen("myColorData.txt", "rb");

	if (pF==NULL) 
	{
		mRet = 1;
		goto mError; 
	}
	else
	{

		do
		{
			mTemp = fgetc(pF);
			printf("%c", mTemp);
			if (mTemp == '~')
				mSize = mFileSize + 1;
			if (mTemp != EOF)
				mFileSize++; 
			
			if (mSize > CHUNK)
			{
				mRet = 4;
				goto mError;
			}
		}while(mTemp!=EOF);
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
	printf("Buffer %s, %p\n", pF, pF);
	return 0;
	mError:
	{
		if (mBuffer != NULL) 
		{
			printf("\nSomething wrong for Free Buffer, try again!\n");
			free(mBuffer);
		}
			
		if (pF != NULL)
		{
			printf("\nSomething wrong for close File, try again!\n");
			fclose(pF);
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

int main(void)
{
	int mRet = 0;
	int mOp = 0;
	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	printf("~                                                            ~\n");
	printf("~                                                            ~\n");
	printf("~  Color Proecss Function                                    ~\n");
	printf("~  Arthor: dgq                                               ~\n");
	printf("~  Version 0.1        Date: Jan 3, 2014                      ~\n");
	printf("~  Hello, listed below is code to process color, enjoy it!!! ~\n");
	printf("~                                                            ~\n");
	printf("~                                                            ~\n");
	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

	//Read the color data from the text file named myColorData.txt
	mRet = mInputDataFromFile();
	if (mRet!=0)
	{
		printf("mInputDataFromFile error!\n");
		return 0;
	}

	do
	{
		printf("\n\n\n0: Test;\n");
		printf("1: Calculate inverse of CCMatrix;\n");
		printf("2: Calculate RGB value before gamma from the value after;\n");
		printf("3: Calculate RGB value after gamma from the value before;\n");
		printf("4: Calculate HSV value from RGB value;\n");
		printf("5: Calculate RGB value from HSV value;\n\n\n");
		printf("\nInput function:  "); 
		scanf("%d", &mOp);
		
		
		switch (mOp) 
		{
		case 0:
			printf("**********************Test E**************************\n");
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
		default:
			printf("***************To be continued loading****************\n");
			break;
		}

	} while ((mOp>=0)&&(mOp<=5));
	return 0;
}

