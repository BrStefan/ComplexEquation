#include<fstream>
#include<cmath>
using namespace std;
class Numar_Complex
{
	private:
		double real,imaginar;
	public: 
		Numar_Complex(double ,double );
		Numar_Complex(const Numar_Complex&);
		~Numar_Complex();
		friend istream& operator >> (istream&, Numar_Complex&);
		friend ostream& operator << (ostream&, const Numar_Complex &);
		friend const Numar_Complex operator + (const Numar_Complex&, const Numar_Complex&);
		friend const Numar_Complex operator - (const Numar_Complex&, const Numar_Complex&);
		friend const Numar_Complex operator * (const Numar_Complex&, const Numar_Complex&);
		friend const Numar_Complex operator / (const Numar_Complex&, const Numar_Complex&);
		friend bool operator > (const Numar_Complex&, const Numar_Complex&);
		friend bool operator == (const Numar_Complex&, const Numar_Complex&);
		friend bool operator != (const Numar_Complex&, const Numar_Complex&);
		void vector(istream&,ostream&);
		void SetReal(double);
		void SetImaginar(double);
		double GetReal();
		double GetImaginar();
		void afisare(ostream&);
		double modul();
		Numar_Complex sqrt();

};

Numar_Complex::Numar_Complex(double r=0, double i=0)
{
	real=r;
	imaginar=i;
}

Numar_Complex::Numar_Complex(const Numar_Complex &nr)
{
	real=nr.real;
	imaginar=nr.imaginar;
}

Numar_Complex::~Numar_Complex()
{
	real=0;
	imaginar=0;
}

istream& operator>> (istream &fin, Numar_Complex &nr)
{
	fin>>nr.real>>nr.imaginar;
	return fin;
}

ostream& operator<< (ostream& fout,const Numar_Complex& nr)
{
	fout<<nr.real<<" "<<nr.imaginar<<"\n";
	return fout;
}
void Numar_Complex::vector(istream &fin , ostream& fout)
{
	fin>>*this;
	fout<<*this;
}
void Numar_Complex::SetReal(double r)
{
	real=r;
}

void Numar_Complex::SetImaginar(double i)
{
	imaginar=i;
}

double Numar_Complex::GetReal()
{
	return real;
}

double Numar_Complex::GetImaginar()
{
	return imaginar;
}

void Numar_Complex::afisare(ostream& fout)
{
	double r,i;

	r=(*this).real;
	i=(*this).imaginar;
	if(r==0 && i==0)fout<<"0";
	else
	{
		if(r)
			fout<<r;
		if(i)
		{
			if(i==1)fout<<"+i";
			else if(i==-1)fout<<"-i";
			else if(i>0)fout<<"+"<<i<<"i";
			else fout<<i<<"i";
		}	
	}
}

double Numar_Complex::modul()
{
	double r,i;
	r=(*this).real;
	i=(*this).imaginar;
	double aux;
	aux=r*r+i*i;
	double modul=pow(aux,0.5);
	return modul;
}

const Numar_Complex operator+(const Numar_Complex& x,const Numar_Complex& y)
{
	Numar_Complex rez;
	rez.imaginar=x.imaginar+y.imaginar;
	rez.real=x.real+y.real;
	return rez;
}

const Numar_Complex operator-(const Numar_Complex& x,const Numar_Complex& y)
{
	Numar_Complex rez;
	rez.imaginar=x.imaginar-y.imaginar;
	rez.real=x.real-y.real;
	return rez;
}

const Numar_Complex operator*(const Numar_Complex& x,const Numar_Complex& y)
{
	Numar_Complex rez;
	rez.real=(x.real*y.real)-(x.imaginar*y.imaginar);
	rez.imaginar=(x.real*y.imaginar)+(x.imaginar*y.real);
	return rez;
}

const Numar_Complex operator/(const Numar_Complex& x,const Numar_Complex& y)
{
	Numar_Complex rez;
	rez.real = ( (x.real*y.real) + (x.imaginar*y.imaginar) ) / ( (y.real*y.real)+(y.imaginar*y.imaginar) );
	rez.imaginar = ( (y.real*x.imaginar) - (y.real*x.imaginar) ) / ( (y.real*y.real)+(y.imaginar*y.imaginar) );
	return rez;
}

bool operator==(const Numar_Complex& x,const Numar_Complex& y)
{
	if(x.real==y.real && x.imaginar==y.imaginar)return true;
	return false;
}

bool operator!=(const Numar_Complex& x,const Numar_Complex& y)
{
	if(x.real==y.real && x.imaginar==y.imaginar)return false;
	return true;
}

bool operator>(const Numar_Complex& x,const Numar_Complex& y)
{
	
	return true;
}

Numar_Complex Numar_Complex::sqrt()
{
	double r,i;
	Numar_Complex rez;

	r=(*this).real;
	i=(*this).imaginar;

	rez.real=pow(r*r+i*i,0.5);
	rez.imaginar=pow(r*r+i*i,0.5);


	rez.real+=r;
	rez.imaginar-=r;

	rez.real/=2;
	rez.imaginar/=2;

	rez.real=pow(rez.real,0.5);
	rez.imaginar=pow(rez.imaginar,0.5);
	return rez;
}
int main()
{
	ifstream fin("date.in");
	ofstream fout("date.out");
	int n,nr,care,care2,care3;
	Numar_Complex *v;
	double cat;
	fin>>n;
	v=new Numar_Complex[n+1];
	fout<<"Compenentele numerelor complexe sunt:\n\n";
	for(int i=1;i<=n;i++)v[i].vector(fin,fout);
		fout<<"\n------------------------\n";
		fin>>nr;
		while(nr)
		{
			switch(nr)
			{
				case 1: // afisare real
				{
					fin>>care;
					if(care<1 || care>n)fout<<"Pozitie invalida\n";
					else
						fout<<"Partea reala a numarului "<<care<<" este: "<<v[care].GetReal()<<"\n";
					break;
				}
				case 2: //afisare imaginar
				{
					fin>>care;
					if(care<1 || care>n)fout<<"Pozitie invalida\n";
					else
						fout<<"Partea imaginara a numarului "<<care<<" este: "<<v[care].GetImaginar()<<"\n";
					break;
				}
				case 3: // seteaza real
				{
					fin>>care>>cat;
					if(care<1 || care>n)fout<<"Pozitie invalida\n";
					else
					{
						v[care].SetReal(cat);
						fout<<"Noua parte reala a numarului "<<care<<" este: "<<cat<<"\n";
					}
					break;
				}
				case 4: // seteaza imaginar
				{
					fin>>care>>cat;
					if(care<1 || care>n)fout<<"Pozitie invalida\n";
					else
					{
						v[care].SetImaginar(cat);
						fout<<"Noua parte imaginara a numarului "<<care<<" este: "<<cat<<"\n";
					}
					break;
				}
				case 5: // afisare modul
				{
					fin>>care;
					if(care<1 || care>n)fout<<"Pozitie invalida\n";
					else
						fout<<"Modulul numarului "<<care<<" este: "<<v[care].modul()<<"\n";
					break;
				}
				case 10: //  afiseaza obiectele
				{
					fout<<"Numerele complexe din lista sunt:\n";
					for(int i=1;i<=n;i++)
						{
							v[i].afisare(fout);
							fout<<"\n";
						}
					break;
				}
				case 6: //  afiseaza radical
				{
					fin>>care;
					if(care<1 || care>n)fout<<"Pozitie invalida\n";
					else
					{	fout<<"Radicalul de ordin 2 al numarului ";
						fout<<care;
						fout<<" este: ";
						Numar_Complex aux;
						aux=v[care].sqrt();
						fout<<aux;
						fout<<"\n";
					}
					break;
				}
				case 7: // adunare
				{
					fin>>care>>care2;
					if(care<1 || care>n)fout<<"Pozitie invalida\n";
					else if(care2<1 || care2>n)fout<<"Pozitie invalida\n";
					else
					{
						Numar_Complex aux;
						aux=v[care]+v[care2];
						fout<<"Numerele ";
						fout<<care;
						fout<<" si ";
						fout<<care2;
						fout<<" adunate dau: ";
						fout<<aux;
						fout<<"\n";
					}
					break;
				}
				case 8: // inmultire
				{
					fin>>care>>care2;
					if(care<1 || care>n)fout<<"Pozitie invalida\n";
					else if(care2<1 || care2>n)fout<<"Pozitie invalida\n";
					else
					{
						Numar_Complex aux;
						aux=v[care]*v[care2];
						fout<<"Numerele ";
						fout<<care;
						fout<<" si ";
						fout<<care2;
						fout<<" inmultite dau: ";
						fout<<aux;
						fout<<"\n";
					}
					break;
				}
				case 9: // impartire
				{
					fin>>care>>care2;
					if(care<1 || care>n)fout<<"Pozitie invalida\n";
					else if(care2<1 || care2>n)fout<<"Pozitie invalida\n";
					else
					{
						if(v[care2]==0)fout<<"Imposibil\n";
						else
						{
							Numar_Complex aux;
							aux=v[care]/v[care2];
							fout<<"Numerele ";
							fout<<care;
							fout<<" si ";
							fout<<care2;
							fout<<" impartite dau: ";
							fout<<aux;
							fout<<"\n";
						}
					}
					break;
				}
				case 11: // ecuatia
				{
					fin>>care>>care2>>care3;
					Numar_Complex aux=0,rez1=0,rez2=0;
					if((care<1 || care>n) || (care2<1 || care2>n) || (care3<1 || care3>n))fout<<"Pozitie invalida\n";
					else if(v[care]==aux && v[care2]!=aux && v[care3]!=aux)
					{
						rez1=v[care3]*(-1);
						rez1=rez1/v[care2];
						fout<<"Rezultatul ecuatiei (";
						v[care2].afisare(fout);
						fout<<")x+(";
						v[care].afisare(fout);
						fout<<") este: "<<rez1;
						fout<<"\n";
					}
					else if(v[care]==aux && v[care2]==aux && v[care3]==aux) fout<<"Invalid\n";
					else
					{
						Numar_Complex delta,r1,r2,delta_rad;
						delta=(v[care2]*v[care2])-4*v[care]*v[care3];
						delta_rad=delta.sqrt();
						r1=((-1)*v[care2]+delta_rad)/2*v[care];
						r2=((-1)*v[care2]-delta_rad)/2*v[care];
						fout<<"Solutia ecuatiei (";
						v[care].afisare(fout);
						fout<<")x^2+(";
						v[care2].afisare(fout);
						fout<<")x+";
						v[care3].afisare(fout);
						fout<<" sunt: ";
						r1.afisare(fout);
						fout<<" si ";
						r2.afisare(fout);
					}
					break;
				}
				default:
				{
					fout<<"Comanda invalida!\n";
					break;
				}
			}
			fin>>nr;
		}

	return 0;
}