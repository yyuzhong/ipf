#include <iostream>
#include <vector>
#include <sys/stat.h>

using namespace std;

long GetFileSize(char* filename)
{
    struct stat stat_buf;
    int rc = stat(filename, &stat_buf);
    return rc == 0 ? stat_buf.st_size : -1;
}


int main(int argc, const char** argv)
{

        char *dataFileName = (char *)argv[1];
        char *labelTrainFile = (char *)argv[2]; 

        long sampleSize = GetFileSize(labelTrainFile)/4;

        printf("%s:%d\n",__FUNCTION__,__LINE__);

        vector<float> trainingLabel(sampleSize);
        vector<float> sampleData(sampleSize);

        FILE *flabel = fopen(labelTrainFile,"r");
        int labelCount=0;
        while( labelCount<sampleSize && (!feof(flabel)) )
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

            //printf("%s:%d, %d:%f\n",__FUNCTION__,__LINE__,labelCount,v);
            trainingLabel[labelCount]= v;
            labelCount++;
        }
        fclose(flabel);

        printf("%s:%d\n",__FUNCTION__,__LINE__);

	FILE *fsample = fopen(dataFileName,"r");
	int sampleCount=0;
	while( sampleCount<sampleSize && (!feof(fsample)) )
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
		printf("%f vs: %f\n",v,trainingLabel[sampleCount]);
		sampleData[sampleCount] = v;
		sampleCount++;
	}
	fclose(fsample);

        printf("%s:%d\n",__FUNCTION__,__LINE__);

        return 0;
}







