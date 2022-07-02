
static double Flag = -1; 
static double k1_rate;
static double k2_rate;

void read_constant(string fname, double& k1)
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
        k1 = stod(line);
        // cout << "c" << i << ": " << line << "\t" << k1 << endl;
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
	read_constant("./k1", k1_rate);
  read_constant("./k2", k2_rate);
  Flag = 1; 

}

double k1Function(double *Value,
                         map <string,int>& NumTrans,
                         map <string,int>& NumPlaces,
                         const vector<string> & NameTrans,
                         const struct InfTr* Trans,
                         const int T,
                         const double& time)
{

    // Definition of the function exploited to calculate the rate,
    // in this case for semplicity we define it throught the Mass Action law
 
 	if( Flag == -1)   init_data_structures();
 
    double intensity = 1.0;
    
    for (unsigned int k=0; k<Trans[T].InPlaces.size(); k++)
    {
      intensity *= pow(Value[Trans[T].InPlaces[k].Id],Trans[T].InPlaces[k].Card);
    }
    
    double rate = (k1_rate/2) * intensity;
    
    return(rate);
}

double k2Function(double *Value,
                  map <string,int>& NumTrans,
                  map <string,int>& NumPlaces,
                  const vector<string> & NameTrans,
                  const struct InfTr* Trans,
                  const int T,
                  const double& time)
  
{
  
  // Definition of the function exploited to calculate the rate,
  // in this case for simplicity we define it through the Mass Action law
  
  if( Flag == -1)   init_data_structures();
  
  double intensity = 1.0;
  
  for (unsigned int k=0; k<Trans[T].InPlaces.size(); k++)
  {
    intensity *= pow(Value[Trans[T].InPlaces[k].Id],Trans[T].InPlaces[k].Card);
  }
  
  double rate = (k2_rate/6) * intensity;
  
  return(rate);
}
