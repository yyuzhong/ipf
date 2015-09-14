// baibo.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
//#include "malloc.h"
//#include "EEJCB11.cpp"
#include <string.h>

//排列个数
#define AN 103
//每个排列上的道数
#define CN 102
//每个道上的采样点个数
#define NN 101
//时窗内采样点个数（奇数）
#define nL 7 //31
#define WS 3


void myread(FILE *flabel,float *bd,int n)
{
    int i = 0;
    while(i<n)
    {
        float v;
        unsigned char temp[4];
        fread((void*)(&temp[0]), sizeof(temp[0]), 1, flabel);
        fread((void*)(&temp[1]), sizeof(temp[1]), 1, flabel);
        fread((void*)(&temp[2]), sizeof(temp[2]), 1, flabel);
        fread((void*)(&temp[3]), sizeof(temp[3]), 1, flabel);
        //int itemp = (temp[0]<<0) | (temp[1]<<8) | (temp[2]<<16) | (temp[3]<<24);
        unsigned int itemp =  (temp[3]<<0) | (temp[2]<<8) | (temp[1]<<16) | (temp[0]<<24);
        v = *((float*)&itemp);
        bd[i] = v;
        i++;
    }
    
}

void B2L(float *bd,int n)
{
/*
    for(int i=0;i<n;i++)
    {
        float elem = bd[i];
        char *temp = (char *)(&elem);
        unsigned int itemp =  (temp[3]<<0) | (temp[2]<<8) | (temp[1]<<16) | (temp[0]<<24);
        float v = *((float*)&itemp);
        bd[i] = v; //did not write back, only for test.
    }
*/
}

void L2B(float *ld){
    /*
    char *temp = (char *)(ld);
    unsigned int itemp =  (temp[3]<<0) | (temp[2]<<8) | (temp[1]<<16) | (temp[0]<<24);
    memcpy((char *)ld,(char *)(&itemp),sizeof(float));
    */
}


//求实对称矩阵的特征值及特征向量的雅格比法
//利用雅格比(Jacobi)方法求实对称矩阵的全部特征值及特征向量
//返回值小于0表示超过迭代jt次仍未达到精度要求
//返回值大于0表示正常返回
//a-长度为n*n的数组，存放实对称矩阵，返回时对角线存放n个特征值
//n-矩阵的阶数
//v-长度为n*n的数组，返回特征向量(按列存储)
//eps-控制精度要求
//jt-整型变量，控制最大迭代次数
int eejcb(float a[],int n,float v[],float eps,int jt)
{ 
	int i,j,p,q,u,w,t,s,l;
	double fm,cn,sn,omega,x,y,d;
	l=1;
	for (i=0; i<=n-1; i++){ 
		v[i*n+i]=1.0;
		for (j=0; j<=n-1; j++){
			if (i!=j) {
				v[i*n+j]=0.0;
			}
		}
	}
	while (1){ 
		fm=0.0;
		for (i=0; i<=n-1; i++){
			for (j=0; j<=n-1; j++){ 
				d=fabs(a[i*n+j]);
				if ((i!=j)&&(d>fm)){ 
					fm=d; 
					p=i; 
					q=j;
				}
			}
		}
		if (fm<eps){
			return(1);
		}
		if (l>jt){
			return(-1);
		}
		l=l+1;
		u=p*n+q; 
		w=p*n+p; 
		t=q*n+p; 
		s=q*n+q;
		x=-a[u];
		y=(a[s]-a[w])/2.0;
		omega=x/sqrt(x*x+y*y);
		if (y<0.0){
			omega=-omega;
		}
		sn=1.0+sqrt(1.0-omega*omega);
		sn=omega/sqrt(2.0*sn);
		cn=sqrt(1.0-sn*sn);
		fm=a[w];
		a[w]=fm*cn*cn+a[s]*sn*sn+a[u]*omega;
		a[s]=fm*sn*sn+a[s]*cn*cn-a[u]*omega;
		a[u]=0.0;
		a[t]=0.0;
		for (j=0; j<=n-1; j++){
			if ((j!=p)&&(j!=q)){ 
				u=p*n+j;
				w=q*n+j;
				fm=a[u];
				a[u]=fm*cn+a[w]*sn;
				a[w]=-fm*sn+a[w]*cn;
			}
		}
		for (i=0; i<=n-1; i++){
			if ((i!=p)&&(i!=q)){ 
				u=i*n+p; 
				w=i*n+q;
				fm=a[u];
				a[u]=fm*cn+a[w]*sn;
				a[w]=-fm*sn+a[w]*cn;
			}
		}
		for (i=0; i<=n-1; i++){ 
			u=i*n+p; 
			w=i*n+q;
			fm=v[u];
			v[u]=fm*cn+v[w]*sn;
			v[w]=-fm*sn+v[w]*cn;
		}
	}
	return(1);
}

//------------------------------------C1基于相关的相干算法------------------------------------------//
//做两道相关
float ccor(float x[],int m,float h[],int n)
{


	float compare(float d[],int m);
	float *y;
	int  nn=m+n-1;
	float Cmax;
	y=(float *)calloc(nn,sizeof(float));
	int i,j;

	for(i=-n+1;i<=m-1;i++)
	{

		y[i+n-1]=0.0;
		for(j=0;j<=m-1;j++)
		{

			if(j-i>=0&&j-i<=n-1)
			{


				y[i+n-1]=y[n+i-1]+x[j]*h[j-i];


			}

		}

	}
	Cmax=compare(y,nn);

	return(Cmax);
	free(y);

}

float compare(float d[],int m)
{

	int i;
	float z;
	z=d[0]; 
	for(i=1;i<m;i++)
		if(d[i]>=z) z=d[i];
	return(z);

}

//做三道相关

float Initial(float a[],float b[],float c[],int n)
{

	float ccor(float x[],int m,float h[],int n);
	float C11,C22,C33,C12,C13;
	float D;// 相干值
	C11=ccor(a,n,a,n);

	C22=ccor(b,n,b,n);

	C33=ccor(c,n,c,n);

	C12=ccor(a,n,b,n);

	C13=ccor(a,n,c,n);

	if((C11==0.0)||(C22==0.0)||(C33==0.0)) 
	{

		D=0.0;

	}

	else D=(float)sqrt((C12/(sqrt(C11*C22)))*(C13/(sqrt(C11*C33))));

        //printf("Init: %f,%f,%f,%f,%f=>%f\n",C11,C22,C33,C12,C13,D);
	return(D);


}

//C1算法
void C1(char* infile, char* outfile)
{


	FILE *fp,*fp1;
	fp=fopen(infile,"rb");
	fp1=fopen(outfile,"wb");
	int i,j,p,q,ii,n;//循环控制变量
	//int k=0;
	/*  float *D1,*D2,*D3;
	    D1=(float *)calloc(N,sizeof(float));
	    D2=(float *)calloc(N,sizeof(float));
	    D3=(float *)calloc(N,sizeof(float));
	    float *a,*b,*c;
	    a=(float *)calloc(nL,sizeof(float));
	    b=(float *)calloc(nL,sizeof(float));
	    c=(float *)calloc(nL,sizeof(float));
	 */float D1[NN], D2[NN], D3[NN];//三道组合
	float a[nL], b[nL], c[nL];//由于计算的时窗内点
	float D1_max,D2_max,D3_max;
	float z;

	for(i=0;i<AN-1;i++)
	{

		for(j=0;j<CN-1;j++)

		{


			fseek(fp,(i*NN*CN+j*NN)*sizeof(float),0);
			//fread(&D1,sizeof(float),NN,fp);
            myread(fp,D1,NN);
			//fread(&D2,sizeof(float),NN,fp);
            myread(fp,D2,NN);
			fseek(fp,(CN-2)*NN*sizeof(float),1);
			//fread(&D3,sizeof(float),NN,fp);
            myread(fp,D3,NN);
			//归一化输入数据                                 
			D1_max=compare(D1,NN);
			D2_max=compare(D2,NN);
			D3_max=compare(D3,NN);
			for(n=0;n<NN;n++)
			{

				D1[n]=D1[n]/D1_max;
				D2[n]=D2[n]/D2_max;
				D3[n]=D3[n]/D3_max;


			}

			for(p=0;p<NN-nL+1;p++)
			{

				for(q=0;q<nL;q++)
				{

					a[q]=D1[p+q];
					b[q]=D2[p+q];
					c[q]=D3[p+q];



				}

				z=Initial(a,b,c,nL);
			        L2B(&z);

				fwrite(&z,sizeof(float),1,fp1);
				//   fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
				// k++;
				//处理时间方向上边界问题
				if(p==0)
				{

					for(ii=0;ii<nL/2;ii++)
					{
						fwrite(&z,sizeof(float),1,fp1);
						//    fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
						//   k++;

					}

				}
				else if(p==NN-nL)
				{

					for(ii=0;ii<nL/2;ii++)
					{
						fwrite(&z,sizeof(float),1,fp1);
						//   fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
						//     k++;


					}

				}
				else continue;


			}
			//处理in-line边界问题 
			if(j==CN-2)
			{


				fseek(fp,(i*NN*CN+j*NN)*sizeof(float),0);
				//fread(&D1,sizeof(float),NN,fp);
                myread(fp,D1,NN);
				//fread(&D2,sizeof(float),NN,fp);
                myread(fp,D2,NN);
				fseek(fp,(CN-2)*NN*sizeof(float),1);
				//fread(&D3,sizeof(float),NN,fp);
                myread(fp,D3,NN);
				//归一化输入数据                                 
				D1_max=compare(D1,NN);
				D2_max=compare(D2,NN);
				D3_max=compare(D3,NN);
				for(n=0;n<NN;n++)
				{

					D1[n]=D1[n]/D1_max;
					D2[n]=D2[n]/D2_max;
					D3[n]=D3[n]/D3_max;


				}

				for(p=0;p<NN-nL+1;p++)
				{

					for(q=0;q<nL;q++)
					{

						a[q]=D1[p+q];



						b[q]=D2[p+q];



						c[q]=D3[p+q];



					}

					z=Initial(a,b,c,nL);
					L2B(&z);

					fwrite(&z,sizeof(float),1,fp1);
					//    fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
					//  k++;
					if(p==0)
					{

						for(ii=0;ii<nL/2;ii++)
						{
							fwrite(&z,sizeof(float),1,fp1);
							//  fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
							//  k++;

						}

					}
					else if(p==NN-nL)
					{

						for(ii=0;ii<nL/2;ii++)
						{
							fwrite(&z,sizeof(float),1,fp1);
							// fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
							// k++;

						}

					}
					else continue;


				}

			}
			else continue;

		}
		//处理cross-line 边界问题
		if(i==AN-2)
		{

			for(j=0;j<CN-1;j++)

			{


				fseek(fp,(i*NN*CN+j*NN)*sizeof(float),0);
				//fread(&D1,sizeof(float),NN,fp);
                myread(fp,D1,NN);
				//fread(&D2,sizeof(float),NN,fp);
                myread(fp,D2,NN);
				fseek(fp,(CN-2)*NN*sizeof(float),1);
				//fread(&D3,sizeof(float),NN,fp);
                myread(fp,D3,NN);
				//归一化输入数据                                 
				D1_max=compare(D1,NN);
				D2_max=compare(D2,NN);
				D3_max=compare(D3,NN);
				for(n=0;n<NN;n++)
				{

					D1[n]=D1[n]/D1_max;
					D2[n]=D2[n]/D2_max;
					D3[n]=D3[n]/D3_max;


				}
				for(p=0;p<NN-nL+1;p++)
				{

					for(q=0;q<nL;q++)
					{

						a[q]=D1[p+q];
						b[q]=D2[p+q];
						c[q]=D3[p+q];

					}

					z=Initial(a,b,c,nL);
					L2B(&z);

					fwrite(&z,sizeof(float),1,fp1);
					//  fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
					// k++;
					if(p==0)
					{

						for(ii=0;ii<nL/2;ii++)
						{
							fwrite(&z,sizeof(float),1,fp1);
							// fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
							//k++;

						}

					}
					else if(p==NN-nL)
					{

						for(ii=0;ii<nL/2;ii++)
						{
							fwrite(&z,sizeof(float),1,fp1);
							//fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
							//   k++;

						}

					}
					else continue;


				}

				if(j==CN-2)
				{


					fseek(fp,(i*NN*CN+j*NN)*sizeof(float),0);
					//fread(&D1,sizeof(float),NN,fp);
                    myread(fp,D1,NN);
					//fread(&D2,sizeof(float),NN,fp);
                    myread(fp,D2,NN);
					fseek(fp,(CN-2)*NN*sizeof(float),1);
					//fread(&D3,sizeof(float),NN,fp);
                    myread(fp,D3,NN);
					//归一化输入数据                                 
					D1_max=compare(D1,NN);
					D2_max=compare(D2,NN);
					D3_max=compare(D3,NN);
					for(n=0;n<NN;n++)
					{

						D1[n]=D1[n]/D1_max;
						D2[n]=D2[n]/D2_max;
						D3[n]=D3[n]/D3_max;


					}
					for(p=0;p<NN-nL+1;p++)
					{

						for(q=0;q<nL;q++)
						{

							a[q]=D1[p+q];



							b[q]=D2[p+q];



							c[q]=D3[p+q];



						}

						z=Initial(a,b,c,nL);
					        L2B(&z);

						fwrite(&z,sizeof(float),1,fp1);
						// fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
						// k++;
						if(p==0)
						{

							for(ii=0;ii<nL/2;ii++)
							{
								fwrite(&z,sizeof(float),1,fp1);
								//  fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
								//  k++;

							}

						}
						else if(p==NN-nL)
						{

							for(ii=0;ii<nL/2;ii++)
							{
								fwrite(&z,sizeof(float),1,fp1);
								//fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
								// k++;

							}

						}
						else continue;


					}

				}
				else continue;

			}


		}
		else continue;

	}
	//free(D1);
	//free(D2);
	//free(D3);
	//free(a);
	//free(b);
	//free(c);

	fclose(fp);
	fclose(fp1);
	//fclose(fp2);

}


//--------------------------------------C2基于相似性的相干算法----------------------------------------//


float Semblance(float a[],float b[],float c[],int n)
{

	float d[nL][WS];//数据矩阵
	float C[WS][WS];//协方差矩阵
	float C_sum=0.0;//C中各元素之和
	float d_sum;//协方差矩阵中的值
	float Tr_C=0.0;     //C的迹
	float D;//相干值
	// printf("%f\n",d[i][0]); printf("%f\n",d[i][0]);
	int i,j,k;
	for(i=0;i<n;i++)
	{

		d[i][0]=a[i];
		d[i][1]=b[i];
		d[i][2]=c[i];


	}
	//求协方差矩阵
	for(i=0;i<WS;i++)
		for(j=0;j<WS;j++)
		{

			d_sum=0.0;
			for(k=0;k<n;k++)
			{


				d_sum+=d[k][i]*d[k][j];



			}
			C[i][j]=d_sum;



		}

	//求C中各元素的和
	for(i=0;i<WS;i++)
		for(j=0;j<WS;j++)
		{

			C_sum+=C[i][j];


		}
	//求矩阵C的迹
	for(i=0;i<WS;i++)
		Tr_C+=C[i][i];

	//求相干值
	if(Tr_C==0.0)
	{

		D=0.0;

	}
	else D=C_sum/(WS*Tr_C);

        //printf("Sem: %f,%f =>%f\n",C_sum,Tr_C,D);

	return(D);


}


void C2(char* infile, char* outfile)
{


	//float Semblance(float a[],float b[],float c[],int nL);
	float compare(float d[],int m);
	FILE *fp,*fp1;
	fp=fopen(infile,"rb");
	fp1=fopen(outfile,"wb");



	int i,j,p,q,ii,n;//循环控制变量
	int k=0;
	float D1[NN],D2[NN],D3[NN];//组合所选的三道
	float a[nL], b[nL], c[nL];//由于计算的时窗内点
	float z;
	float D1_max,D2_max,D3_max;
	for(i=0;i<AN-1;i++)
	{

		for(j=0;j<CN-1;j++)

		{


			fseek(fp,(i*NN*CN+j*NN)*sizeof(float),0);
			//fread(&D1,sizeof(float),NN,fp);
            myread(fp,D1,NN);
			//fread(&D2,sizeof(float),NN,fp);
            myread(fp,D2,NN);
			fseek(fp,(CN-2)*NN*sizeof(float),1);
			//fread(&D3,sizeof(float),NN,fp);
            myread(fp,D3,NN);
			//归一化输入数据                                 
			D1_max=compare(D1,NN);
			D2_max=compare(D2,NN);
			D3_max=compare(D3,NN);
            //printf("Semblance: %d: %f,%f,%f\n",__LINE__,D1_max,D2_max,D3_max);
			for(n=0;n<NN;n++)
			{

				D1[n]=D1[n]/D1_max;
				D2[n]=D2[n]/D2_max;
				D3[n]=D3[n]/D3_max;


			}
			for(p=0;p<NN-nL+1;p++)
			{

				for(q=0;q<nL;q++)
				{

					a[q]=D1[p+q];
					b[q]=D2[p+q];
					c[q]=D3[p+q];



				}




				z=Semblance(a,b,c,nL);
			        L2B(&z);


				fwrite(&z,sizeof(float),1,fp1);
				//     fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
				// k++;
				//处理时间方向上边界问题
				if(p==0)
				{

					for(ii=0;ii<nL/2;ii++)
					{
						fwrite(&z,sizeof(float),1,fp1);
						//          fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
						//            k++;

					}

				}
				else if(p==NN-nL)
				{

					for(ii=0;ii<nL/2;ii++)
					{
						fwrite(&z,sizeof(float),1,fp1);
						//        fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
						//              k++;


					}

				}
				else continue;


			}
			//处理in-line边界问题 
			if(j==CN-2)
			{

				fseek(fp,(i*NN*CN+j*NN)*sizeof(float),0);
				//fread(&D1,sizeof(float),NN,fp);
                myread(fp,D1,NN);
				//fread(&D2,sizeof(float),NN,fp);
                myread(fp,D2,NN);
				fseek(fp,(CN-2)*NN*sizeof(float),1);
				//fread(&D3,sizeof(float),NN,fp);
                myread(fp,D3,NN);
				//归一化输入数据                                 
				D1_max=compare(D1,NN);
				D2_max=compare(D2,NN);
				D3_max=compare(D3,NN);
                //printf("Semblance: %d: %f,%f,%f\n",__LINE__,D1_max,D2_max,D3_max);
				for(n=0;n<NN;n++)
				{

					D1[n]=D1[n]/D1_max;
					D2[n]=D2[n]/D2_max;
					D3[n]=D3[n]/D3_max;


				}
				for(p=0;p<NN-nL+1;p++)
				{

					for(q=0;q<nL;q++)
					{

						a[q]=D1[p+q];



						b[q]=D2[p+q];



						c[q]=D3[p+q];



					}





					z=Semblance(a,b,c,nL);
			                L2B(&z);

					fwrite(&z,sizeof(float),1,fp1);
					//     fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
					//  k++;
					if(p==0)
					{

						for(ii=0;ii<nL/2;ii++)
						{
							fwrite(&z,sizeof(float),1,fp1);
							//  fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
							//  k++;

						}

					}
					else if(p==NN-nL)
					{

						for(ii=0;ii<nL/2;ii++)
						{
							fwrite(&z,sizeof(float),1,fp1);
							// fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
							//          k++;

						}

					}
					else continue;


				}

			}
			else continue;

		}
		//处理cross-line 边界问题
		if(i==AN-2)
		{

			for(j=0;j<CN-1;j++)

			{


				fseek(fp,(i*NN*CN+j*NN)*sizeof(float),0);
				//fread(&D1,sizeof(float),NN,fp);
                myread(fp,D1,NN);
				//fread(&D2,sizeof(float),NN,fp);
                myread(fp,D2,NN);
				fseek(fp,(CN-2)*NN*sizeof(float),1);
				//fread(&D3,sizeof(float),NN,fp);
                myread(fp,D3,NN);
				//归一化输入数据                                 
				D1_max=compare(D1,NN);
				D2_max=compare(D2,NN);
				D3_max=compare(D3,NN);
                //printf("Semblance: %d: %f,%f,%f\n",__LINE__,D1_max,D2_max,D3_max);
				for(n=0;n<NN;n++)
				{

					D1[n]=D1[n]/D1_max;
					D2[n]=D2[n]/D2_max;
					D3[n]=D3[n]/D3_max;


				}
				for(p=0;p<NN-nL+1;p++)
				{

					for(q=0;q<nL;q++)
					{

						a[q]=D1[p+q];
						b[q]=D2[p+q];
						c[q]=D3[p+q];

					}


					z=Semblance(a,b,c,nL);
			                L2B(&z);

					fwrite(&z,sizeof(float),1,fp1);
					//    fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
					// k++;
					if(p==0)
					{

						for(ii=0;ii<nL/2;ii++)
						{
							fwrite(&z,sizeof(float),1,fp1);
							//   fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
							//           k++;

						}

					}
					else if(p==NN-nL)
					{

						for(ii=0;ii<nL/2;ii++)
						{
							fwrite(&z,sizeof(float),1,fp1);
							//    fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
							//             k++;

						}

					}
					else continue;


				}

				if(j==CN-2)
				{


					fseek(fp,(i*NN*CN+j*NN)*sizeof(float),0);
					//fread(&D1,sizeof(float),NN,fp);
                    myread(fp,D1,NN);
					//fread(&D2,sizeof(float),NN,fp);
                    myread(fp,D2,NN);
					fseek(fp,(CN-2)*NN*sizeof(float),1);
					//fread(&D3,sizeof(float),NN,fp);
                    myread(fp,D3,NN);
					//归一化输入数据                                 
					D1_max=compare(D1,NN);
					D2_max=compare(D2,NN);
					D3_max=compare(D3,NN);

                    //printf("Semblance: %d: %f,%f,%f\n",__LINE__,D1_max,D2_max,D3_max);

					for(n=0;n<NN;n++)
					{

						D1[n]=D1[n]/D1_max;
						D2[n]=D2[n]/D2_max;
						D3[n]=D3[n]/D3_max;

					}

					for(p=0;p<NN-nL+1;p++)
					{

						for(q=0;q<nL;q++)
						{

							a[q]=D1[p+q];



							b[q]=D2[p+q];



							c[q]=D3[p+q];



						}


						z=Semblance(a,b,c,nL);
			                        L2B(&z);

						fwrite(&z,sizeof(float),1,fp1);
						//     fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
						// k++;
						if(p==0)
						{

							for(ii=0;ii<nL/2;ii++)
							{
								fwrite(&z,sizeof(float),1,fp1);
								//  fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
								//        k++;

							}

						}
						else if(p==NN-nL)
						{

							for(ii=0;ii<nL/2;ii++)
							{
								fwrite(&z,sizeof(float),1,fp1);
								//    fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
								//            k++;

							}

						}
						else continue;


					}

				}
				else continue;

			}


		}
		else continue;  

	}
	//free(D1);
	//free(D2);
	//free(D3);
	//free(a);
	//free(b);
	//free(c);

	fclose(fp);
	fclose(fp1);
	//fclose(fp2);

}

//---------------------------------C3基于特征值的相干算法------------------------------------//


//求矩阵的最大特征值
float Eigenvalue(float a[],float b[],float c[],int n)
{

	float compare(float d[],int m);
	//int eejcb(float *a,int n,float *v,float eps,int jt);

	float v[WS][WS];
	float Eigen_v[WS];//存放所得特征值
	float E_max;//特征值最大值
	float eps;

	eps=0.001f;
	int ii;
	float  d[nL][WS];//数据矩阵
	float C[WS][WS];//协方差矩阵
	float C_sum=0.0;//C中各元素之和
	float  d_sum;//协方差矩阵中的值
	float Tr_C=0.0;     //C的迹
	float D;//相干值

	int i,j,k;
	for(i=0;i<n;i++)
	{

		d[i][0]=a[i];
		d[i][1]=b[i];
		d[i][2]=c[i];


	}
	//求协方差矩阵
	for(i=0;i<WS;i++)
		for(j=0;j<WS;j++)
		{

			d_sum=0.0;
			for(k=0;k<n;k++)
			{


				d_sum+=d[k][i]*d[k][j];



			}
			C[i][j]=d_sum;

		}
	ii=eejcb(C[0],WS,v[0],eps,100);

	if (ii>0)
	{

		for (ii=0; ii<WS; ii++)
			Eigen_v[ii]=C[ii][ii];


	}
	E_max=compare(Eigen_v,WS);


	//求矩阵C的迹
	for(i=0;i<WS;i++)
		Tr_C+=C[i][i];

	//求相干值
	if(Tr_C==0.0)
	{

		D=0.0;

	}
	else D=E_max/Tr_C;

	return(D);


}




void C3(char* infile, char* outfile)
{


	FILE *fp,*fp1;
	fp=fopen(infile,"rb");
	fp1=fopen(outfile,"wb");

	int i,j,p,q,ii,n;//循环控制变量
	//int k=0;
	/*  float *D1,*D2,*D3;
	    D1=(float *)calloc(N,sizeof(float));
	    D2=(float *)calloc(N,sizeof(float));
	    D3=(float *)calloc(N,sizeof(float));
	    float *a,*b,*c;
	    a=(float *)calloc(nL,sizeof(float));
	    b=(float *)calloc(nL,sizeof(float));
	    c=(float *)calloc(nL,sizeof(float));
	 */float D1[NN], D2[NN], D3[NN];//三道组合
	float a[nL], b[nL], c[nL];//由于计算的时窗内点
	float D1_max,D2_max,D3_max;
	float z;

	for(i=0;i<AN-1;i++)
	{

		for(j=0;j<CN-1;j++)

		{


			fseek(fp,(i*NN*CN+j*NN)*sizeof(float),0);
			//fread(&D1,sizeof(float),NN,fp);
            myread(fp,D1,NN);
			//fread(&D2,sizeof(float),NN,fp);
            myread(fp,D2,NN);
			fseek(fp,(CN-2)*NN*sizeof(float),1);
			//fread(&D3,sizeof(float),NN,fp);
            myread(fp,D3,NN);
			//归一化输入数据                                 
			D1_max=compare(D1,NN);
			D2_max=compare(D2,NN);
			D3_max=compare(D3,NN);
			for(n=0;n<NN;n++)
			{

				D1[n]=D1[n]/D1_max;
				D2[n]=D2[n]/D2_max;
				D3[n]=D3[n]/D3_max;


			}
			for(p=0;p<NN-nL+1;p++)
			{

				for(q=0;q<nL;q++)
				{

					a[q]=D1[p+q];
					b[q]=D2[p+q];
					c[q]=D3[p+q];



				}

				z=Eigenvalue(a,b,c,nL);
			        L2B(&z);

				fwrite(&z,sizeof(float),1,fp1);
				//     fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
				// k++;
				//处理时间方向上边界问题
				if(p==0)
				{

					for(ii=0;ii<nL/2;ii++)
					{
						fwrite(&z,sizeof(float),1,fp1);
						//           fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
						//             k++;

					}

				}
				else if(p==NN-nL)
				{

					for(ii=0;ii<nL/2;ii++)
					{
						fwrite(&z,sizeof(float),1,fp1);
						//        fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
						//              k++;


					}

				}
				else continue;


			}
			//处理in-line边界问题 
			if(j==CN-2)
			{


				fseek(fp,(i*NN*CN+j*NN)*sizeof(float),0);
				//fread(&D1,sizeof(float),NN,fp);
                myread(fp,D1,NN);
				//fread(&D2,sizeof(float),NN,fp);
                myread(fp,D2,NN);
				fseek(fp,(CN-2)*NN*sizeof(float),1);
				//fread(&D3,sizeof(float),NN,fp);
                myread(fp,D3,NN);
				//归一化输入数据                                 
				D1_max=compare(D1,NN);
				D2_max=compare(D2,NN);
				D3_max=compare(D3,NN);
				for(n=0;n<NN;n++)
				{

					D1[n]=D1[n]/D1_max;
					D2[n]=D2[n]/D2_max;
					D3[n]=D3[n]/D3_max;


				}
				for(p=0;p<NN-nL+1;p++)
				{

					for(q=0;q<nL;q++)
					{

						a[q]=D1[p+q];



						b[q]=D2[p+q];



						c[q]=D3[p+q];



					}

					z=Eigenvalue(a,b,c,nL);
			                L2B(&z);

					fwrite(&z,sizeof(float),1,fp1);
					//      fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
					//  k++;
					if(p==0)
					{

						for(ii=0;ii<nL/2;ii++)
						{
							fwrite(&z,sizeof(float),1,fp1);
							//  fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
							//  k++;

						}

					}
					else if(p==NN-nL)
					{

						for(ii=0;ii<nL/2;ii++)
						{
							fwrite(&z,sizeof(float),1,fp1);
							//   fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
							//            k++;

						}

					}
					else continue;


				}

			}
			else continue;

		}
		//处理cross-line 边界问题
		if(i==AN-2)
		{

			for(j=0;j<CN-1;j++)

			{


				fseek(fp,(i*NN*CN+j*NN)*sizeof(float),0);
				//fread(&D1,sizeof(float),NN,fp);
                myread(fp,D1,NN);
				//fread(&D2,sizeof(float),NN,fp);
                myread(fp,D2,NN);
				fseek(fp,(CN-2)*NN*sizeof(float),1);
				//fread(&D3,sizeof(float),NN,fp);
                myread(fp,D3,NN);
				//归一化输入数据                                 
				D1_max=compare(D1,NN);
				D2_max=compare(D2,NN);
				D3_max=compare(D3,NN);
				for(n=0;n<NN;n++)
				{

					D1[n]=D1[n]/D1_max;
					D2[n]=D2[n]/D2_max;
					D3[n]=D3[n]/D3_max;


				}
				for(p=0;p<NN-nL+1;p++)
				{

					for(q=0;q<nL;q++)
					{

						a[q]=D1[p+q];
						b[q]=D2[p+q];
						c[q]=D3[p+q];

					}

					z=Eigenvalue(a,b,c,nL);
					L2B(&z);

					fwrite(&z,sizeof(float),1,fp1);
					//    fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
					// k++;
					if(p==0)
					{

						for(ii=0;ii<nL/2;ii++)
						{
							fwrite(&z,sizeof(float),1,fp1);
							//   fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
							//           k++;

						}

					}
					else if(p==NN-nL)
					{

						for(ii=0;ii<nL/2;ii++)
						{
							fwrite(&z,sizeof(float),1,fp1);
							//    fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
							//            k++;

						}

					}
					else continue;


				}

				if(j==CN-2)
				{


					fseek(fp,(i*NN*CN+j*NN)*sizeof(float),0);
					//fread(&D1,sizeof(float),NN,fp);
                    myread(fp,D1,NN);
					//fread(&D2,sizeof(float),NN,fp);
                    myread(fp,D2,NN);
					fseek(fp,(CN-2)*NN*sizeof(float),1);
					//fread(&D3,sizeof(float),NN,fp);
                    myread(fp,D3,NN);
					//归一化输入数据                                 
					D1_max=compare(D1,NN);
					D2_max=compare(D2,NN);
					D3_max=compare(D3,NN);
					for(n=0;n<NN;n++)
					{

						D1[n]=D1[n]/D1_max;
						D2[n]=D2[n]/D2_max;
						D3[n]=D3[n]/D3_max;


					}
					for(p=0;p<NN-nL+1;p++)
					{

						for(q=0;q<nL;q++)
						{

							a[q]=D1[p+q];



							b[q]=D2[p+q];



							c[q]=D3[p+q];



						}

						z=Eigenvalue(a,b,c,nL);
					        L2B(&z);

						fwrite(&z,sizeof(float),1,fp1);
						//     fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
						// k++;
						if(p==0)
						{

							for(ii=0;ii<nL/2;ii++)
							{
								fwrite(&z,sizeof(float),1,fp1);
								// fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
								//        k++;

							}

						}
						else if(p==NN-nL)
						{

							for(ii=0;ii<nL/2;ii++)
							{
								fwrite(&z,sizeof(float),1,fp1);
								//    fprintf(fp2,"%d, %d, %d, %d ,%f\n",i,j,p,k,z);
								//           k++;

							}

						}
						else continue;


					}

				}
				else continue;

			}


		}
		else continue;  

	}
	//free(D1);
	//free(D2);
	//free(D3);
	//free(a);
	//free(b);
	//free(c);

	fclose(fp);
	fclose(fp1);


}



int main()
{
    
	int ID;//选择算法
	/* char *string;
     string = "There are five coherency algorithms to caculate coherence cube.They are respectively:\n";
     printf("%s\n",string);
     string = "C1--Initial Coherency Algorithm;\n";
     printf("%s\n",string);
     string = "C2--Semblance-based Coherency Algorithm;\n";
     printf("%s\n",string);
     string = "C3--EigenStructure-based Coherency Algorithm;\n";
     printf("%s\n",string);
     string = "C4--Local Structural Entropy Coherency Algorithm;\n";
     printf("%s\n",string);
     string = "C5--Higher-Order-Statistics-based Coherence-estimation method;\n";
     printf("%s\n",string);
     
     string = "These five algorithms are below.Please choose the algorithm by which you want\nto calculate the coherence cube.\n";
     printf("%s\n",string);
     
     string = "1-----------------C1\n";
     printf("%s\n",string);
     string = "2-----------------C2\n";
     printf("%s\n",string);
     string = "3-----------------C3\n";
     printf("%s\n",string);
     string = "4-----------------C4\n";
     printf("%s\n",string);
     string = "5-----------------C5\n";
     printf("%s\n",string);
	 */ 

     printf("%s\n", "Input your choise\n");
     scanf("%d", &ID);
    
    if(ID==1) C1("/data/seis/fah/gx.dat","/data/seis/fah/ch1.dat");
    if(ID==2) C2("/data/seis/fah/gx.dat","/data/seis/fah/ch2.dat");
    if(ID==3) C3("/data/seis/fah/gx.dat","/data/seis/fah/ch3.dat");
    
    return(ID);
}




/*

   }
   D=compare(z,nn);
   return (D);
   free(z);
   free(bb_corr);
   free(cc_corr);
   free(abc_corr);



   }
 */
