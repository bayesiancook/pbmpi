#include "Partition.h"
#include "BiologicalSequences.h"
#include "StringStreamUtils.h"
#include <algorithm>
#include <map>

using namespace std;

vector<PartitionScheme> PartitionProcess::ReadSchemes(string schemefile, int Nsite, int myid, bool linkgam, bool unlinkgtr, string rrtype)
{
	string error = "Error: improperly formatted scheme file\n";

	ifstream theStream(schemefile.c_str());
	if (! theStream)   {
		cerr << "error: cannot open partition scheme file: " << schemefile << "\n";
		exit(1);
	}

	map<string, size_t> fixrrparts;
	map<string, size_t> fixstatparts;

	PartitionScheme rrscheme(Nsite);
	PartitionScheme statscheme(Nsite);
	PartitionScheme dgamscheme(Nsite);

	vector<size_t> partsites;

	// get the unlinked partition info
	string line;
	size_t lineNum = 0;
	while(getline(theStream, line))
	{
	    if(line.find_first_not_of(" \t\n\v\f\r") == std::string::npos)
	    {
	        continue;
	    }

	    if(line.find_last_of(";") != std::string::npos)
	        line.erase(line.find_last_of(";"));

		stringstream ss(line);

		string item;
		vector<string> elems;
		while (getline(ss, item, ',')) {
			elems.push_back(item);
		}

		if(elems.size() == 1)
		{
			cerr << error;
			cerr << "Line " << lineNum << "\n";
			exit(1);
		}
		else if(elems.size() == 0)
		{
		    continue;
		}

		string type;
		bool fixprof;
		bool estimate = false;

		for(vector<string>::iterator it = elems.begin(); it != elems.end(); it++)
		{
			if(it == elems.begin())
			{
				ss.clear();
				ss.str(*it);

				ss >> type;

				if(type.empty())
				{
						cerr << error;
						cerr << "Line " << lineNum << ": partition type not specified\n";
						exit(1);
				}

				std::transform(type.begin(), type.end(), type.begin(), ::tolower);

				fixprof = ((*(type.rbegin()) != 'f' && *(type.rbegin()) != 'x') || type == "dayhoff" );
				estimate = *(type.rbegin()) == 'x';

				if(!fixprof)
				{
					// remove 'f' or 'e' character
					if (type.size () > 0)  type.resize (type.size () - 1);
				}

				if(type == "gtr")
				{
					type = "None";
					if(unlinkgtr)
					{
						stringstream ss(type);
						ss << rrscheme.Npart;
						type = ss.str();
					}

					if(fixrrparts.find(type) == fixrrparts.end())
					{
						rrscheme.partType.push_back("None");
						fixrrparts[type] = rrscheme.Npart++;
					}

					if(fixprof)
					{
						fixprof = false;
						estimate = true;
					}
				}
				else
				{
					if(fixprof && fixstatparts.find(type) == fixstatparts.end())
					{
						statscheme.partType.push_back(type);
						fixstatparts[type] = statscheme.Npart++;
					}

					if(fixrrparts.find(type) == fixrrparts.end())
					{
						rrscheme.partType.push_back(type);
						fixrrparts[type] = rrscheme.Npart++;
					}
				}
			}
			else
			{
				ss.clear();
				ss.str(*it);

				while(ss >> item) {}

				ss.clear();
				ss.str(item);

				vector<int> range;
				while (getline(ss, item, '-'))
				{
					if(IsInt(item))
					{
						range.push_back(atoi(item.c_str()) - 1);
					}
					else
					{
						cerr << error;
						cerr << "Line " << lineNum << ": token '" << item << "'\n";
						exit(1);
					}
				}

				if(range.empty() || range.size() > 2)
				{
					cerr << error;
					cerr << "Line " << lineNum << "\n";
					exit(1);
				}

				if(range.size() == 1)
				{
					partsites.push_back(range[0]);

					dgamscheme.sitePart[range[0]] = dgamscheme.Npart;

					rrscheme.sitePart[range[0]] = fixrrparts[type];

					if(fixprof)
					{
						statscheme.sitePart[range[0]] = fixstatparts[type];
					}
					else
					{
						statscheme.sitePart[range[0]] = statscheme.Npart;
					}
				}
				else
				{
					for(size_t i = range[0]; i <= range[1]; i++)
					{
						partsites.push_back(i);

						dgamscheme.sitePart[i] = dgamscheme.Npart;

						rrscheme.sitePart[i] = fixrrparts[type];

						if(fixprof)
						{
							statscheme.sitePart[i] = fixstatparts[type];
						}
						else
						{
							statscheme.sitePart[i] = statscheme.Npart;
						}
					}
				}
			}
		}

		dgamscheme.Npart++;

		if(!fixprof)
		{
			if(estimate)
				statscheme.partType.push_back("None");
			else
				statscheme.partType.push_back("Empirical");
			statscheme.Npart++;
		}

		lineNum++;
	}

	if(partsites.size() != Nsite)
	{
		size_t rrpart;
		if(fixrrparts.find("None") != fixrrparts.end())
		{
			rrpart = fixrrparts["None"];
		}
		else
		{
			rrscheme.partType.push_back("None");
			rrpart = rrscheme.Npart++;
		}

		statscheme.partType.push_back("None");

		std::sort(partsites.begin(), partsites.end());

		vector<size_t>::iterator it = partsites.begin();
		for(size_t site = 0; site < Nsite; site++)
		{
			if(*it != site)
			{
				rrscheme.sitePart[site] = rrpart;
				statscheme.sitePart[site] = statscheme.Npart;
				dgamscheme.sitePart[site] = dgamscheme.Npart;
			}
			else
			{
				it++;
			}
		}

		statscheme.Npart++;
		dgamscheme.Npart++;
	}

	if(linkgam)
	{
		dgamscheme.Npart = 1;
		std::fill(dgamscheme.sitePart.begin(), dgamscheme.sitePart.end(), 0);
	}

	if(rrtype != "Part")
	{
		if(rrtype == "None" && unlinkgtr)
		{
			for(size_t p = 0; p < rrscheme.partType.size(); p++)
			{
				rrscheme.partType[p] = "None";
			}
		}
		else
		{
			rrscheme.Npart = 1;
			std::fill(rrscheme.sitePart.begin(), rrscheme.sitePart.end(), 0);
			rrscheme.partType.clear();
			rrscheme.partType.push_back(rrtype);
		}
	}

	rrscheme.update();
	statscheme.update();
	dgamscheme.update();

	vector<PartitionScheme> schemes;

	schemes.push_back(rrscheme);
	schemes.push_back(statscheme);
	schemes.push_back(dgamscheme);

	if(myid == 0)
	{
		stringstream ss;
		ss << "Partition scheme\n";
		map<string, size_t> stat_type_counts;
		map<string, size_t> rr_type_counts;

		size_t max_label = 0;
		size_t max_count = 0;
		size_t offset = 20;
		for(size_t i = 0; i < schemes[0].Npart; i++)
		{
			string t = schemes[0].partType[i] == "None" ? "gtr" : schemes[0].partType[i];

			rr_type_counts[t]++;

			max_label = std::max(t.size(), max_label);
			max_count = std::max(rr_type_counts[t], max_count);
		}
		for(size_t i = 0; i < schemes[1].Npart; i++)
		{
			string t = schemes[1].partType[i] == "None" ? "Free" : schemes[1].partType[i];

			stat_type_counts[t]++;

			max_label = std::max(t.size(), max_label);
			max_count = std::max(stat_type_counts[t], max_count);
		}

		string s = "Relative rates: ";
		for(map<string,size_t>::iterator it = rr_type_counts.begin(); it != rr_type_counts.end(); it++)
		{
			ss << right << setw(offset) << s << left << setw(max_label + 1) << it->first << setw(max_count + 1) << it->second << " partitions\n";
			if(it == rr_type_counts.begin())
				s = "";
		}
		s = "Profiles: ";
		for(map<string,size_t>::iterator it = stat_type_counts.begin(); it != stat_type_counts.end(); it++)
		{
			ss << right << setw(offset) << s << left << setw(max_label + 1) << it->first << setw(max_count + 1) << it->second << " partitions\n";
			if(it == stat_type_counts.begin())
				s = "";
		}

		ss << setw(offset) << "Rates across sites:";
		ss << schemes[2].Npart << " partitions" << endl << endl;
		cerr << ss.str();
	}

	return schemes;
}

