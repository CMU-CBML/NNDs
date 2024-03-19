#include "BasicDataStructure.h"

Edge::Edge()
{
	act=1;
	pt[0]=-1; pt[1]=-1;
	//face[0]=-1; face[1]=-1;
	len=1.;
	//pned[0]=-1; pned[1]=-1;
	//pnfc[0]=-1; pnfc[1]=-1;
	pn[0][0]=3; pn[1][0]=3;//0 for edge, 1 for face (1st edge diretion), 2 for face (2nd edge direction), 3 for end
	pn[0][1]=-1; pn[1][1]=-1;//edge or face ID
	prt=-1;
	chd[0]=-1; chd[0]=-1;
	midpt=-1;

	add = -1;
	dual = 0;
	//pdual[0] = 0.; pdual[1] = 0.; pdual[2] = 0.;
	ploc = -1;
	spoke = 0;
}

bool Edge::operator==(const Edge& ed)
{
	return ((pt[0]==ed.pt[0] && pt[1]==ed.pt[1])||(pt[0]==ed.pt[1] && pt[1]==ed.pt[0]));
}

/////////////////////////////////////////////////////

Vertex::Vertex()
{
	coor[0]=0.;	coor[1]=0.;	coor[2]=0.;
	coortmp[0]=0.;	coortmp[1]=0.;	coortmp[2]=0.;
	w=1.;
	wtmp=0.;
	update=0;
	//ctmp[0]=0.; ctmp[1]=0.; ctmp[2]=0.;
	type=0;
	trun=0;
	truntmp=0;
	act=1;
	//cal=0;
	//lev=0;
	//h2a=-1;
	//sc=1.;
	//id=0;
	//vln=4;
	//ntr=16;
	//vn.clear();
	//sup.clear();
	//child.clear();
	//c.clear();
	//r1.clear();
	//r2.clear();
	//err=0.;
	kitvU[0]=1.; kitvU[1]=1.; kitvU[2]=1.; kitvU[3]=1.;
	kitvV[0]=1.; kitvV[1]=1.; kitvV[2]=1.; kitvV[3]=1.;
}

Vertex Vertex::operator+(const Vertex& v)
{
	Vertex a;
	a.coor[0]=coor[0]+v.coor[0];
	a.coor[1]=coor[1]+v.coor[1];
	a.coor[2]=coor[2]+v.coor[2];
	return a;
}

Vertex Vertex::operator-(const Vertex& v)
{
	Vertex a;
	a.coor[0]=coor[0]-v.coor[0];
	a.coor[1]=coor[1]-v.coor[1];
	a.coor[2]=coor[2]-v.coor[2];
	return a;
}

Vertex Vertex::operator*(double a)
{
	Vertex v;
	v.coor[0]=coor[0]*a;
	v.coor[1]=coor[1]*a;
	v.coor[2]=coor[2]*a;
	return v;
}

Vertex Vertex::operator/(double a)
{
	Vertex v;
	v.coor[0]=coor[0]/a;
	v.coor[1]=coor[1]/a;
	v.coor[2]=coor[2]/a;
	return v;
}

bool Vertex::operator==(const Vertex& v)
{
	double tol(0.00001);
	double dif=sqrt((coor[0]-v.coor[0])*(coor[0]-v.coor[0])+(coor[1]-v.coor[1])*(coor[1]-v.coor[1])+(coor[2]-v.coor[2])*(coor[2]-v.coor[2]));
	//return (lev==v.lev&&coor[0]==v.coor[0]&&coor[1]==v.coor[1]&&coor[2]==v.coor[2]);
	return (lev==v.lev && dif<tol);
}

////////////////////////////////////////////////////////////

Element::Element()
{
	cnct[0]=-1; cnct[1]=-1; cnct[2]=-1; cnct[3]=-1;
	edge[0]=-1; edge[1]=-1; edge[2]=-1; edge[3]=-1;
	//edge.resize(4);
	act=0;
	//aff=0;
	lev=0;
	lv=0.;
	h2a=-1;
	type=0;
	//square=1;
	//nxp=0;
	prt=-1;
	chd[0]=-1; chd[1]=-1; chd[2]=-1; chd[3]=-1;
	//trun=0;
	//err=0.;
	ref=0;
	nT=0;
	Tjunc[0]=-1; Tjunc[1]=-1; Tjunc[2]=-1; Tjunc[3]=-1;
	//nbed[0]=-1; nbed[1]=-1; nbed[2]=-1; nbed[3]=-1;
	pn[0][0]=3; pn[1][0]=3; pn[2][0]=3; pn[3][0]=3;//0 for edge, 1 for face (1st edge diretion), 2 for face (2nd edge direction), 3 for end
	pn[0][1]=-1; pn[1][1]=-1; pn[2][1]=-1; pn[3][1]=-1;
	//edge_nb[0][0]=-1; edge_nb[1][0]=-1; edge_nb[2][0]=-1; edge_nb[3][0]=-1;
	//edge_nb[0][1]=-1; edge_nb[1][1]=-1; edge_nb[2][1]=-1; edge_nb[3][1]=-1;
	//corn_nb[0]=-1; corn_nb[1]=-1; corn_nb[2]=-1; corn_nb[3]=-1;
	newpID.clear();
	//newpID[0]=-1; newpID[1]=-1; newpID[2]=-1; newpID[3]=-1; newpID[4]=-1;

	//parent[0]=-1;
	//parent[1]=-1;
	//edge.resize(4);
	//IEN.resize(16);
	//UV.resize(16);
	kv.clear();
	IEN.clear();
	//edge_neibors.resize(4);
	//corner_neibors.resize(4);
	parent[0]=-1; parent[1]=-1;
	child.clear();

	add = -1;
	dual = 0;
	focus = 0;
	//pdual[0] = 0.; pdual[1] = 0.; pdual[2] = 0.;
	cnctd[0] = 0; cnctd[1] = 0; cnctd[2] = 0; cnctd[3] = 0;
	//edged[0] = 0; edged[1] = 0; edged[2] = 0; edged[3] = 0;
	edlen[0] = 1.; edlen[1] = 1.; edlen[2] = 1.; edlen[3] = 1.;

	update = 0;
	c1 = 0;
}

bool Element::operator==(const Element& e)
{
	if(lev!=e.lev /*|| vertx_neibors.size()!=e.vertx_neibors.size()*/)
		return false;
	else
	{
		vector<int> copy1(cnct,cnct+4),copy2(e.cnct,e.cnct+4);
		sort(copy1.begin(),copy1.end());
		sort(copy2.begin(),copy2.end());
		for(unsigned int i=0;i<copy1.size();i++)
		{
			if(copy1[i]!=copy2[i])
			{
				return false;
			}
		}
		return true;
	}
}

bool PointIPP::operator==(const PointIPP& pt)
{
	if(index[0]==pt.index[0] && index[1]==pt.index[1] && pm[0]==pt.pm[0] && pm[1]==pt.pm[1])
	{
		return true;
	}
	else
	{
		return false;
	}
}

////////////////////////////////////////////////////////////////

BezierElement::BezierElement()
{
	order=3;
	for(int i=0;i<16;i++)
	{
		pts[i][0]=0.; pts[i][1]=0.; pts[i][2]=0.;
		w[i]=0.;
	}
	for(int i=0;i<25;i++)
	{
		pts4[i][0]=0.; pts4[i][1]=0.; pts4[i][2]=0.;
		w4[i]=0.;
	}
	type = 0;
	focus = 0;

	for (int i = 0; i < 4; i++)
	{
		bc[i] = 0;
	}
	prt = -1;
}

void BezierElement::Basis(double u, double v, double Nt[16], double dNdt[16][2]) const
{
	double Nu[4]={(1.-u)*(1.-u)*(1.-u),3.*(1.-u)*(1.-u)*u,3.*(1.-u)*u*u,u*u*u};
	double Nv[4]={(1.-v)*(1.-v)*(1.-v),3.*(1.-v)*(1.-v)*v,3.*(1.-v)*v*v,v*v*v};
	double dNdu[4]={-3.*(1.-u)*(1.-u),3.-12.*u+9.*u*u,3.*(2.-3.*u)*u,3.*u*u};
	double dNdv[4]={-3.*(1.-v)*(1.-v),3.-12.*v+9.*v*v,3.*(2.-3.*v)*v,3.*v*v};
	int i,j,loc(0);
	for(i=0;i<4;i++)
	{
		for(j=0;j<4;j++)
		{
			Nt[loc]=Nu[j]*Nv[i];
			dNdt[loc][0]=dNdu[j]*Nv[i];
			dNdt[loc][1]=Nu[j]*dNdv[i];
			loc++;
		}
	}
}

void BezierElement::Basis4(double u, double v, double Nt[25], double dNdt[25][2]) const
{
	double Nu[5]={(1.-u)*(1.-u)*(1.-u)*(1.-u),4.*(1.-u)*(1.-u)*(1.-u)*u,6.*(1.-u)*(1.-u)*u*u,4.*(1.-u)*u*u*u,u*u*u*u};
	double Nv[5]={(1.-v)*(1.-v)*(1.-v)*(1.-v),4.*(1.-v)*(1.-v)*(1.-v)*v,6.*(1.-v)*(1.-v)*v*v,4.*(1.-v)*v*v*v,v*v*v*v};
	double dNdu[5]={-4.*(1.-u)*(1.-u)*(1.-u),4.*(1.-u)*(1.-u)*(1.-4.*u),12.*u*(1.-3.*u+2.*u*u),4.*(3.-4.*u)*u*u,4.*u*u*u};
	double dNdv[5]={-4.*(1.-v)*(1.-v)*(1.-v),4.*(1.-v)*(1.-v)*(1.-4.*v),12.*v*(1.-3.*v+2.*v*v),4.*(3.-4.*v)*v*v,4.*v*v*v};
	int i,j,loc(0);
	for(i=0;i<5;i++)
	{
		for(j=0;j<5;j++)
		{
			Nt[loc]=Nu[j]*Nv[i];
			dNdt[loc][0]=dNdu[j]*Nv[i];
			dNdt[loc][1]=Nu[j]*dNdv[i];
			loc++;
		}
	}
}

void BezierElement::Para2Phys(double u, double v, double pt[3])
{
	double Nx[16];
	double dNdt[16][2];
	Basis(u,v,Nx,dNdt);
	pt[0]=0.; pt[1]=0.; pt[2]=0.;
	for(int i=0;i<16;i++)
	{
		pt[0]+=pts[i][0]*Nx[i];
		pt[1]+=pts[i][1]*Nx[i];
		pt[2]+=pts[i][2]*Nx[i];
	}
}

void BezierElement::Para2Phys4(double u, double v, double pt[3])
{
	double Nx[25];
	double dNdt[25][2];
	Basis4(u,v,Nx,dNdt);
	pt[0]=0.; pt[1]=0.; pt[2]=0.;
	for(int i=0;i<25;i++)
	{
		pt[0]+=pts4[i][0]*Nx[i];
		pt[1]+=pts4[i][1]*Nx[i];
		pt[2]+=pts4[i][2]*Nx[i];
	}
}

void BezierElement::SurfPointNormal(double u, double v, array<double,3>& pt, array<double,3>& nm) const
{
	double Nx[16];
	double dNdt[16][2];
	Basis(u,v,Nx,dNdt);
	pt[0]=0.; pt[1]=0.; pt[2]=0.;
	nm[0]=0.; nm[1]=0.; nm[2]=0.;
	double nmtmp[2][3]={{0.,0.,0.},{0.,0.,0.}};
	for(int i=0;i<16;i++)
	{
		pt[0]+=pts[i][0]*Nx[i];
		pt[1]+=pts[i][1]*Nx[i];
		pt[2]+=pts[i][2]*Nx[i];
		nmtmp[0][0]+=pts[i][0]*dNdt[i][0];
		nmtmp[0][1]+=pts[i][1]*dNdt[i][0];
		nmtmp[0][2]+=pts[i][2]*dNdt[i][0];
		nmtmp[1][0]+=pts[i][0]*dNdt[i][1];
		nmtmp[1][1]+=pts[i][1]*dNdt[i][1];
		nmtmp[1][2]+=pts[i][2]*dNdt[i][1];
	}
	nm[0]=nmtmp[0][1]*nmtmp[1][2]-nmtmp[0][2]*nmtmp[1][1];
	nm[1]=nmtmp[0][2]*nmtmp[1][0]-nmtmp[0][0]*nmtmp[1][2];
	nm[2]=nmtmp[0][0]*nmtmp[1][1]-nmtmp[0][1]*nmtmp[1][0];
	double len=sqrt(nm[0]*nm[0]+nm[1]*nm[1]+nm[2]*nm[2]);
	nm[0]/=len; nm[1]/=len; nm[2]/=len;
}

void BezierElement::SurfPointNormal4(double u, double v, array<double,3>& pt, array<double,3>& nm) const
{
	double Nx[25];
	double dNdt[25][2];
	Basis4(u,v,Nx,dNdt);
	pt[0]=0.; pt[1]=0.; pt[2]=0.;
	nm[0]=0.; nm[1]=0.; nm[2]=0.;
	double nmtmp[2][3]={{0.,0.,0.},{0.,0.,0.}};
	for(int i=0;i<25;i++)
	{
		pt[0]+=pts4[i][0]*Nx[i];
		pt[1]+=pts4[i][1]*Nx[i];
		pt[2]+=pts4[i][2]*Nx[i];
		nmtmp[0][0]+=pts4[i][0]*dNdt[i][0];
		nmtmp[0][1]+=pts4[i][1]*dNdt[i][0];
		nmtmp[0][2]+=pts4[i][2]*dNdt[i][0];
		nmtmp[1][0]+=pts4[i][0]*dNdt[i][1];
		nmtmp[1][1]+=pts4[i][1]*dNdt[i][1];
		nmtmp[1][2]+=pts4[i][2]*dNdt[i][1];
	}
	nm[0]=nmtmp[0][1]*nmtmp[1][2]-nmtmp[0][2]*nmtmp[1][1];
	nm[1]=nmtmp[0][2]*nmtmp[1][0]-nmtmp[0][0]*nmtmp[1][2];
	nm[2]=nmtmp[0][0]*nmtmp[1][1]-nmtmp[0][1]*nmtmp[1][0];
	double len=sqrt(nm[0]*nm[0]+nm[1]*nm[1]+nm[2]*nm[2]);
	nm[0]/=len; nm[1]/=len; nm[2]/=len;
}

/////////////////////////

void RegularPatchBasis::Evaluate(double x,double a,double b)//input a(t-b)
{
	double t=a*(x-b);
	if(t<0.||t>1.)
	{
		for(int i=0;i<4;i++) {val[i]=0.; Dval[i]=0.;}
	}
	else
	{
		val[0]=(1.-3.*t+3.*t*t-t*t*t)/6.;
		val[1]=(4.-6.*t*t+3.*t*t*t)/6.;
		val[2]=(1.+3.*t+3.*t*t-3.*t*t*t)/6.;
		val[3]=t*t*t/6.;
		Dval[0]=a*(-3.+6.*t-3.*t*t)/6.;
		Dval[1]=a*(-12.*t+9.*t*t)/6.;
		Dval[2]=a*(3.+6.*t-9.*t*t)/6.;
		Dval[3]=a*t*t/2.;
		D2val[0]=a*(6.-6.*t)/6.;
		D2val[1]=a*(-12.+18.*t)/6.;
		D2val[2]=a*(6.-18.*t)/6.;
		D2val[3]=a*2.*t/2.;
	}
}

void RegularPatchBasis::Clear()
{
	for(int i=0;i<4;i++)
	{
		val[i]=0.;
		Dval[i]=0.;
		D2val[i]=0.;
	}
}