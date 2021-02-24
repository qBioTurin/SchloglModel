
static double Flag = -1; 
static double p_c1;
static double p_c2;
static double p_c3;
static double p_c4;

void read_constant(string fname, double& p)
{
	ifstream f (fname);
	string line;
	if(f.is_open())
	{
		int i = 1;
		while (getline(f,line))
		{
			switch(i)
			{
				case 1:
					p = stod(line);
					
					break;
			}
			++i;
		}
		f.close();
    }
	else
	{
		std::cerr<<"\nUnable to open " << fname << ": file do not exists\": file do not exists\n";
		exit(EXIT_FAILURE);
	}
}

void init_data_structures()
{
	read_constant("./c1", p_c1);
    read_constant("./c2", p_c2);
    read_constant("./c3", p_c3);
    read_constant("./c4", p_c4);
    cout << "p_c1: " << p_c1 << "\n" << endl;
    cout << "p_c2: " << p_c2 << "\n" << endl;
    cout << "p_c3: " << p_c3 << "\n" << endl;
    cout << "p_c4: " << p_c4 << "\n" << endl;
    Flag = 1; 

}

double c1(double *Value,
                         map <string,int>& NumTrans,
                         map <string,int>& NumPlaces,
                         const vector<string> & NameTrans,
                         const struct InfTr* Trans,
                         const int T,
                         const double& time)
{

    // Definition of the function exploited to calculate the rate,
    // in this case for semplicity we define it throught the Mass Action  law
 
 	if( Flag == -1)   init_data_structures();
 
    double intensity = 1;
    
	for (unsigned int k=0; k<Trans[T].InPlaces.size(); k++)
	{
		intensity *= pow(Value[Trans[T].InPlaces[k].Id],Trans[T].InPlaces[k].Card);
        //cout << "p_c1: " << p_c1 << " intensity: "<< intensity << " place:  " << Value[Trans[T].InPlaces[k].Id]<<"\n" << endl;
	}
	
	double rate = p_c1/2 * intensity;

    return(rate);
}

double c2(double *Value,
                         map <string,int>& NumTrans,
                         map <string,int>& NumPlaces,
                         const vector<string> & NameTrans,
                         const struct InfTr* Trans,
                         const int T,
                         const double& time)
{

    // Definition of the function exploited to calculate the rate,
    // in this case for semplicity we define it throught the Mass Action  law
 
 	if( Flag == -1)   init_data_structures();
 
    double intensity = 1;
    
	for (unsigned int k=0; k<Trans[T].InPlaces.size(); k++)
	{
		intensity *= pow(Value[Trans[T].InPlaces[k].Id],Trans[T].InPlaces[k].Card);
        //cout << "p_c2: " << p_c2 << " intensity: "<< intensity << " place:  " << Value[Trans[T].InPlaces[k].Id]<<"\n" << endl;
	}
	
	double rate = p_c2/6 * intensity;
    
    return(rate);
}

double c3(double *Value,
                         map <string,int>& NumTrans,
                         map <string,int>& NumPlaces,
                         const vector<string> & NameTrans,
                         const struct InfTr* Trans,
                         const int T,
                         const double& time)
{

    // Definition of the function exploited to calculate the rate,
    // in this case for semplicity we define it throught the Mass Action  law
 
 	if( Flag == -1)   init_data_structures();
 
    double intensity = 1;
    
	for (unsigned int k=0; k<Trans[T].InPlaces.size(); k++)
	{
		intensity *= pow(Value[Trans[T].InPlaces[k].Id],Trans[T].InPlaces[k].Card);
	}
	
	double rate = p_c3 * intensity;

    return(rate);
}

double c4(double *Value,
                         map <string,int>& NumTrans,
                         map <string,int>& NumPlaces,
                         const vector<string> & NameTrans,
                         const struct InfTr* Trans,
                         const int T,
                         const double& time)
{

    // Definition of the function exploited to calculate the rate,
    // in this case for semplicity we define it throught the Mass Action  law
 
 	if( Flag == -1)   init_data_structures();
 
    double intensity = 1;
    
	for (unsigned int k=0; k<Trans[T].InPlaces.size(); k++)
	{
		intensity *= pow(Value[Trans[T].InPlaces[k].Id],Trans[T].InPlaces[k].Card);
	}
	
	double rate = p_c4 * intensity;

    return(rate);
}
