#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<cmath>
#include<cstring>
#include<vector>

using namespace std;

struct lammpsdata
{
    int id;
    double r;   //Variable to store the radius in Polar co-ordinates
    double th;    //Variable to store theta in Polar co-ordinates
    double z;   //Variable to store z in Polar co-ordinates
    double vr;  //Variable to store the velocity along r
    double vth;   //Variable to store the velocity along theta
    double vz;    //Variable to store the velocity along z

};

// To store global information provided by a dump file
struct dumpInfo
{
  string filename; // Stores the path to the file
  int first_timestep; // Stores the first timestep 
  int last_timestep; // Stores the last timestep
  int dump_frequency; // Stores the difference b/n adjacent timesteps
  int number_of_atoms; // Stores the number of atoms. This is assumed to be the same in all the timesteps in a dump file. 
  vector<string> atom_info; // Stores the names of all atom information 

  // Stores the coordinates of the bounding box
  double xlo; 
  double xhi;
  double ylo;
  double yhi;
  double zlo;
  double zhi;

  // Stores the nature of boundaries
  string boundX;
  string boundY;
  string boundZ;
};

// Stores data in one timestep
struct timestepData
{
  vector< vector <double> > data;
};

// Function to estimate the size of a file in bytes

// Arguments :
// char* file - (Input) - Path to file

//Return value :
// size of the file in bytes

streampos fileSize(const char* file)
{
    streampos fsize = 0;
    ifstream inp(file, ios::in);

    fsize = inp.tellg();
    inp.seekg(0, ios::end);
    fsize = inp.tellg() - fsize;
    inp.close();
    return fsize;
}

// A general method to read data and store information in a dump file

// Arguments :
// char* file - (Input) - Path to the file
// struct dumpInfo* info - (Output) - the struct where the information will be written

// Return value :
// 1 - if file is read successfully
// 0 - if it encounters some error.

int getDumpInfo(const char* file, struct dumpInfo* info)
{
  bool first_step = true;
  bool second_step = false;
  int step;
  double size;
  
  string line;
  string temp;

  ifstream inp (file);
  
  if(!inp.is_open())
    {
      cout<<"Error opening the file :"<<file<<endl;
      return 0;
    }
  size = fileSize(file);
  info->filename = file;
  while(getline(inp, line))
    {
      istringstream iss(line);
      iss>>temp;
      if(strcmp(temp.c_str(),"ITEM:")==0)
	{
	  iss>>temp;
	  if(strcmp(temp.c_str(),"TIMESTEP")==0)
	    {
	      inp>>step;
	      if(first_step)
		info->first_timestep = step;
	      if(second_step)
		{
		  info->dump_frequency = step - info->first_timestep;
		  second_step = false;
		  inp.clear();
		  inp.seekg(size - 1.3*(info->number_of_atoms * static_cast<int>(info->atom_info.size()) * sizeof(double)), ios::beg); // Multiplying 1.3 here is a cheap trick. We may have to think of something better.
		}
	    }
	  else if(strcmp(temp.c_str(),"NUMBER")==0 && first_step)
	    {
	      iss>>temp>>temp;
	      inp>>info->number_of_atoms;	
	    }
	  else if(strcmp(temp.c_str(),"BOX")==0 && first_step)
	    {
	      iss>>temp;
	      iss>>info->boundX>>info->boundY>>info->boundZ;
	      inp>>info->xlo>>info->xhi;
	      inp>>info->ylo>>info->yhi;
	      inp>>info->zlo>>info->zhi;
	    }
	  else if(strcmp(temp.c_str(),"ATOMS")==0 && first_step)
	    {
	      while(iss>>temp)
		{
		  info->atom_info.push_back(temp);
		}
	      second_step=true;
	      first_step=false;
	      inp.clear();
	      inp.seekg((info->number_of_atoms * static_cast<int>(info->atom_info.size()) * sizeof(double)), ios::beg);
	    }
	}
    }
  info->last_timestep=step;
  cout<<"First Timestep:"<<info->first_timestep<<endl;
  cout<<"Last Timestep:"<<info->last_timestep<<endl;
  cout<<"Dump Frequency:"<<info->dump_frequency<<endl;
  cout<<"X Bounds:"<<info->xlo<<" "<<info->xhi<<endl;
  cout<<"Y Bounds:"<<info->ylo<<" "<<info->yhi<<endl;
  cout<<"Z Bounds:"<<info->zlo<<" "<<info->zhi<<endl;
  cout<<"Number of atoms:"<<info->number_of_atoms<<endl;
  cout<<"Atom info:\t";
  for(int i = 0;i < static_cast<int>(info->atom_info.size());i++)
    cout<<info->atom_info[i]<<" ";
  cout<<endl;
  inp.close();
  return 1;
}

int main()
{
    string readLine;
    int timestep, natoms;
    struct dumpInfo info;
    
    /*Initializing the velocity sums and avereage variables for the entire timesteps*/
    double VrTimeAvg=0.0;
    double VzTimeAvg=0.0;

    /*Getting info about the number of timesteps*/
    cout<<endl<<"Enter the number of timesteps: ";
    cin>>timestep;

    /*Getting the data for the bins*/
    double r1, r2, z1, z2;
    int rdiv, zdiv;
    cout<<endl<<"Enter the limits of r for the binning - r1: ";
    cin>>r1;
    cout<<"\t r2: ";
    cin>>r2;
    cout<<"\n No: of divisions along the limits of r: ";
    cin>>rdiv;
    cout<<endl<<"Enter the limits of z for the binning - z1: ";
    cin>>z1;
    cout<<"\t z2: ";
    cin>>z2;
    cout<<"\n No: of divisions along the limits of z: ";
    cin>>zdiv;

    double** VrTimeSum=new double*[zdiv];
    for(int i=0; i<zdiv; i++)
        VrTimeSum[i]=new double[rdiv];

    double** VzTimeSum=new double*[zdiv];
    for(int j=0; j<zdiv; j++)
        VzTimeSum[j]=new double[rdiv];
    
    /*Finding the lower and upper limits of r and z*/
    double zlow=min(z1,z2);
    double rlow=min(r1,r2);
    double zup=z1+z2-zlow;
    double rup=r1+r2-rlow;
    double zl1=zlow;
    double zu1;
    double rl1=rlow;
    double ru1; //Defining the limits to be used within the loop//
    ifstream inp("dump.flow");  /*Opening the dump file to read the data*/
    ofstream out("vflow.txt");  /*Opening a file to write the velocity for each bins for the single timesteps*/
    
    getDumpInfo("dump.flow", &info);
    
    /*Loop for processing the data from the file */
    for(int z=0; z<timestep; z++)
    {
        for(int p=0; p<7; p++)
        inp>>readLine;
        inp>>natoms;
        cout<<endl<<"No of atoms: \t"<<natoms<<endl;
        for(int p=0; p<21; p++)
            inp>>readLine;
        lammpsdata* ld= new lammpsdata[natoms]; //Allocating structure using 'new'
        double x,y,vx,vy,tht;
    
        /*Storing the data from the Lammps dump file into the data structure*/
      for(int i=0; i<natoms; i++)
        {
            inp>>ld[i].id;
            inp>>x;
            inp>>y;
            inp>>ld[i].z;
            inp>>vx;
            inp>>vy;
            inp>>ld[i].vz;
            ld[i].r=sqrt((x*x)+(y*y));
            ld[i].th=(atan(y/x));
            tht=(atan(y/x));
            ld[i].vr=(vx*cos(tht))+(vy*sin(tht));
            ld[i].vth=-vx*sin(tht)+vy*cos(tht);

        }

        double vrSum=0.0, vzSum=0.0;        //Declaring and initializing the velocity sum variables//
        double vrAvg, vzAvg;
        int counter=0;
        out<<"\n r(lower)\t r(upper)\t z(lower)\t z(upper)\t Vr(avg)\t Vz(avg)";
        zu1=zup;
        zl1=zlow;
        for(int a=0; a<zdiv; a++)       //Loop for z limits division //
        {
            zu1=zl1+((zup-zlow)/zdiv);
            rl1=rlow;
            ru1=rup;
            for(int b=0; b<rdiv; b++)       //Loop for r limits division//
            {
                ru1=rl1+((rup-rlow)/rdiv);
                vrSum=0.0;
                vzSum=0.0;
                counter=0;
                for(int c=0; c<natoms; c++)
                {
                    if(ld[c].r<=ru1 && ld[c].r>=rl1 && ld[c].z<=zu1 && ld[c].z>=zl1)      //Checking whether the atom is within the limits//
                    {   vrSum+=ld[c].vr;
                        vzSum+=ld[c].vz;
                        counter++;

                    }

                }
                if(counter==0)
                {
                    vrAvg=0;
                    vzAvg=0;
                }
                else
                {
                    vrAvg=vrSum/counter;
                    vzAvg=vzSum/counter;
                }
                out<<"\n"<<rl1<<"\t"<<ru1<<"\t"<<zl1<<"\t"<<zu1<<"\t"<<vrAvg<<"\t"<<vzAvg;
                VrTimeSum[a][b]+=vrAvg;
                VzTimeSum[a][b]+=vzAvg;
                rl1=ru1;
            }
            zl1=zu1;
        }
            delete []ld;
        }


    out.close();
    inp.close();
    ofstream outfinal("vTimeAvg.txt");
    zu1=zup;
    zl1=zlow;
    outfinal<<"\n r(lower)\t r(upper)\t z(lower)\t z(upper)\t r(mid)\t z(mid)\t Vr(time avg)\t Vz(time avg)";
    for(int l=0; l<zdiv; l++)
    {
        zu1=zl1+((zup-zlow)/zdiv);
        rl1=rlow;
        ru1=rup;
        for(int m=0; m<rdiv; m++)
        {
            ru1=rl1+((rup-rlow)/rdiv);
            double rMid=(ru1+rl1)/2;
            double zMid=(zu1+zl1)/2;
            VrTimeAvg=VrTimeSum[l][m]/timestep;
            VzTimeAvg=VzTimeSum[l][m]/timestep;
            outfinal<<"\n"<<rl1<<"\t"<<ru1<<"\t"<<zl1<<"\t"<<zu1<<"\t"<<rMid<<"\t"<<zMid<<"\t"<<VrTimeAvg<<"\t"<<VzTimeAvg;
            rl1=ru1;
        }
        zl1=zu1;
    }
    outfinal.close();
    for(int i=0; i<zdiv; i++)
        delete [] VrTimeSum[i];
    delete [] VrTimeSum;
    for(int i=0; i<zdiv; i++)
        delete [] VzTimeSum[i];
    delete [] VzTimeSum;
    
    return 0;
}
