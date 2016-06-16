#include "SequenceAlignment.h"

int main(int argc, char* argv[])	{

	SequenceAlignment* data = new FileSequenceAlignment(argv[1],1,0);
	ofstream os(argv[2]);
	data->ToStream(os);

}
