
#ifndef PROTEINSEQUENCEALIGNMENT_H
#define PROTEINSEQUENCEALIGNMENT_H


class ProteinSequenceAlignment : public SequenceAlignment	{

	public:

	ProteinSequenceAlignment(CodonSequenceAlignment* from)	{
		Ntaxa = from->Ntaxa;
		Nsite = from->Nsite;
		taxset = from->taxset;
		statespace = new ProteinStateSpace();

		Data = new int*[Ntaxa];
		for (int i=0; i<Ntaxa; i++)	{
			Data[i] = new int[Nsite];
			for (int j=0; j<Nsite; j++)	{
				int tmp = from->GetState(i,j);
				if (tmp == unknown)	{
					Data[i][j] = unknown;
				}
				else	{
					Data[i][j] = from->GetCodonStateSpace()->Translation(tmp);
				}
			}
		}
	}
};

#endif

