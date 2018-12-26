
#include "VoronoiCurve.h"

namespace Voronoicrv
{

VD::Face_handle farr[10000];
float bb[10000][5];
int fin=0,bi=0;


bool selected_face(VD::Face_handle fit1)
{
    for(int i=0;i<fin;i++)
        if(fit1==farr[i])
            return true;
    return false;
}

int VoronoiCurve::iterator1(VD::Vertex_iterator vit2)
{
    int i;
    double l11[2],l12[2],l21[2],l22[2];
    DT::Vertex_handle vit3;
    VD::Halfedge_around_vertex_circulator he1=vit2->incident_halfedges();
    VD::Face_handle fh1;
    VD::Halfedge_handle hh1;
    int ub=0;
    for(i=0;i<3;i++)
    {
        if(he1->face()->is_unbounded()||selected_face(he1->face()))
        {
            ub++;
            if(ub==3)
                return 0;
        }
        else
        {
            fh1=he1->face();
            vit3=he1->face()->dual();
            farr[fin]=he1->face();
            fin++;
        }
        he1++;
    }
    float x1,x2,y1,y2;
    x1=vit2->point().x();
    y1=vit2->point().y();
    x2=vit3->point().x();
    y2=vit3->point().y();
    double m=(y2-y1)/(x2-x1);
    if(m!=0)
    {
        double b1=((1/m)*x2)+y2;
        l11[0]=x1;
        l11[1]=y1;
        l12[0]=x2;
        l12[1]=y2;
        l21[0]=x2-10000;
        l22[0]=x2+10000;
        l21[1]=(-(1/m)*(x2-10000))+b1;
        l22[1]=(-(1/m)*(x2+10000))+b1;
        Point_2 p1=Point_2(x2-100,(-(1/m)*(x2-100))+b1);
        Point_2 p2=Point_2(x2+100,(-(1/m)*(x2+100))+b1);
        int lt=0;
        if(CGAL::left_turn(p1,p2,vit2->point()))
            lt=1;
        VD::Ccb_halfedge_circulator chc1=fh1->halfedge()->ccb();
        VD::Ccb_halfedge_circulator chc2=chc1;
        for(i=0;i<1000;i++)
        {
            if(lt==0&&!CGAL::left_turn(p1,p2,chc1->source()->point()))
                vd.set_status(chc1->source()->point());
            else
                if(lt==1&&CGAL::left_turn(p1,p2,chc1->source()->point()))
                {
                    vd.set_status(chc1->source()->point());
                    for(int k1=bi-1;k1>=0;k1--)
                    {
                        bb[k1+1][0]=bb[k1][0];
                        bb[k1+1][1]=bb[k1][1];
                    }
                    bb[0][0]=chc1->source()->point().x();
                    bb[0][1]=chc1->source()->point().y();
                    bi++;
                }
            chc1--;
            if(chc1->source()->point()==chc2->source()->point())
                break;
        }}
    return 1;
}

VD::Vertex_iterator VoronoiCurve::find_it(Point_2 p)
{
    int a,b,c,d;
    VD::Vertex_iterator vit1=vd.vertices_begin();
    do{
        if(abs((int)vit1->point().x()-(int)p.x())<0.00001&&abs((int)vit1->point().y()-(int)p.y())<0.00001)
            return vit1;
        vit1++;
    }while(vit1!=vd.vertices_end());
    printf("Errorrrrrrrrrrr\n in %f %f",p.x(),p.y());
    return NULL;
}


int VoronoiCurve::getIndex(double x, double y)
{
	for(int i=0; i<(int)points.size(); ++i)
	{
		if ((points[i].first==x)&&(points[i].second==y))
			return i;
	}
	
	return -1;		
}

void VoronoiCurve::collectBoundaryIndices(int pci, int cind, Point_2 pcon[][2], Point_2 con[][2])
{
	int numEdges=0;
    //consider only outer peel
        for(int i=0;i<pci;i++)
        {     
             _boundary.push_back(pair<int, int>(getIndex(pcon[i][0].x(),pcon[i][0].y()), getIndex(pcon[i][1].x(),pcon[i][1].y())));				
            numEdges++;
			
            
        }    
   
        for(int i=0;i<cind;i++)
        {
            {
                _boundary.push_back(pair<int, int>(getIndex(con[i][0].x(),con[i][0].y()), getIndex(con[i][1].x(),con[i][1].y())));				
            numEdges++;
            }
        }
		
	_boundary.resize(numEdges);
}



void VoronoiCurve::reconstruct()
{   
int ci=0,li=0,pci=0,mi=0,emi=0,cind=0;
double chull[1000][2];
Point_2 mat[1000][2],pcon[1000][2],emat[1000][2],con[1000][2];



    	Site_2 site;	
	int i, count1;

	for (i = 0; i < (int)points.size(); i++)
	{
		site= Site_2(points[i].first, points[i].second);
		vd.insert(site);
	}
    
    int kl=0;

bloop1:
    i=0;
    Point_2 p;
    int n,unb=0,b=0;
    n=i;
    Face_iterator fit=vd.faces_begin();
    do{
        if(fit->is_unbounded())
            unb++;
        else
            b++;
        fit++;
    }while(fit!=vd.faces_end());
    int bou=0;
    int n1c=0;
    DT::Vertex_circulator  		vc=vd.dual().incident_vertices(vd.dual().infinite_vertex()),done(vc);
    do{
        chull[ci][0]=vc->point().x();
        chull[ci][1]=vc->point().y();
        ci++;
    }while(++vc!=done);
    int allconvex=0;
    DT::Face_iterator fit1=vd.dual().faces_begin();
    if(ci==points.size())
    {
        allconvex=1;
        n1c=1;
        goto complete1;
    }

    
    do{
        int fl=0;
        for(int i=0;i<3;i++){
            if(vd.dual().is_infinite(fit1->neighbor(i)))
                if(fl==0)
                {
                    bou++;
                    fl=1;
                }
        }
        fit1++;
    }while(fit1!=vd.dual().faces_end());

complete1:
    VD::Halfedge_iterator hit=vd.halfedges_begin();
//    if(n1c==0)
    {
    DT::Face_iterator fit2=vd.dual().faces_begin();
    float big=0;
    float farthest[2];
	
	//Extract the extremal Voronoi vertices
    do{
        DT::Face_handle f;
        if (vd.dual().is_infinite(vd.dual().locate((vd.dual().circumcenter(fit2)))))
        {
            float d=abs(sqrt((vd.dual().triangle(fit2).vertex(0).x()-vd.dual().circumcenter(fit2).x())*(vd.dual().triangle(fit2).vertex(0).x()-vd.dual().circumcenter(fit2).x())
                             +(vd.dual().triangle(fit2).vertex(0).y()-vd.dual().circumcenter(fit2).y())*(vd.dual().triangle(fit2).vertex(0).y()-vd.dual().circumcenter(fit2).y())));
            bb[bi][0]=vd.dual().circumcenter(fit2).x();
            bb[bi][1]=vd.dual().circumcenter(fit2).y();
            bb[bi][2]=d;
            bb[bi][3]=0;
            if(d>big)
            {
                d=big;
                farthest[0]=bb[bi][0];
                farthest[1]=bb[bi][1];
            }
            bi++;
        }
        
        fit2++;
    }while(fit2!=vd.dual().faces_end());
	
	
    VD::Vertex_iterator vit1=vd.vertices_begin();
    VD::Vertex_iterator vit2;
    do{
        if(vit1->point().x()==farthest[0]&&vit1->point().y()==farthest[1])
            vit2=vit1;
        vit1++;
    }while(vit1!=vd.vertices_end());
	
	//sorting the extremal Voronoi vertices
    count1=vd.initialize_status();
    for(i=0;i<bi;i++)
        for(int j=i;j<bi;j++)
        {
            if(bb[i][2]<bb[j][2])
            {
                int tarr[4];
                tarr[0]=bb[j][0];
                tarr[1]=bb[j][1];
                tarr[2]=bb[j][2];
                tarr[3]=bb[j][3];
                bb[j][0]=bb[i][0];
                bb[j][1]=bb[i][1];
                bb[j][2]=bb[i][2];
                bb[j][3]=bb[i][3];
                bb[i][0]=tarr[0];
                bb[i][1]=tarr[1];
                bb[i][2]=tarr[2];
                bb[i][3]=tarr[3];
            }
        }
    Point_2 p1;
    i=0;
    while(1)
    {
        p1=Point_2(bb[i][0],bb[i][1]);
    nloop:
        
        if(abs(p1.x())<0.0000001&&abs(p1.y())<0.0000001)
            goto comp;
        int ri=iterator1(find_it(p1));
        vd.set_visit(p1);
        vit1=vd.vertices_begin();
        VD::Halfedge_around_vertex_circulator he2=find_it(p1)->incident_halfedges();
        int fl1=0;
        for(int j=0;j<3;j++)
        {
            if(he2->has_source())
            {
                if(vd.check_visit(he2->source()->point())==0&&vd.check_status(he2->source()->point())==0)
                {
                    if(fl1==1)
                    {
                        bb[bi][0]=he2->source()->point().x();
                        bb[bi][1]=he2->source()->point().y();
                        bi++;
                        goto njmp;
                    }
                    else
                        p1=he2->source()->point();
                    fl1=1;
                }
            }
            he2++;
        }
		
    njmp:
        if(fl1==1)
            goto nloop;
        if(ri==0)
        {
            for(int k=0;k<bi;k++)
                if(vd.check_visit(Point_2(bb[k][0],bb[k][1]))==0)
                {
                    p1=Point_2(bb[k][0],bb[k][1]);
                    if(bb[k][0]==0.000&&bb[k][1]==0.0000)
                        goto comp;
                    goto nloop;
                }
            goto comp;
        }
    }
	
comp:
    VD::Vertex_iterator vit4=vd.vertices_begin();
    for(i=0;i<count1;i++)
        vit4++;
    
    hit=vd.halfedges_begin();
    do{
        if(hit->has_source()&&hit->has_target())
            if((vd.check_status(hit->source()->point())!=vd.check_status(hit->target()->point())))
            {
                pcon[pci][0]=hit->face()->dual()->point();;
                pcon[pci][1]=hit->opposite()->face()->dual()->point();
                pci++;
            }
        hit++;
    }while(hit!=vd.halfedges_end());
	
	
    Point_2 l1,l2;
    hit=vd.halfedges_begin();
    do{
        if(!hit->has_source()&&hit->has_target())
        {
            if(vd.check_status(hit->target()->point()))
            {
                con[cind][0]=hit->face()->dual()->point();
                con[cind][1]=hit->opposite()->face()->dual()->point();
                cind++;
            }
        }
        if(hit->has_source()&&!hit->has_target())
        {
            if(vd.check_status(hit->source()->point()))
            {
                con[cind][0]=hit->face()->dual()->point();
                con[cind][1]=hit->opposite()->face()->dual()->point();
                cind++;
            }
        }
        if(((!hit->has_source()&&hit->has_target()))||((hit->has_source()&&!hit->has_target())))
        {
            if((hit->has_source()))
            {
                if(vd.check_status(hit->source()->point()))
                {
                    l1=hit->face()->dual()->point();
                    l2=hit->opposite()->face()->dual()->point();
                    int fl=1;
                    for(int yi=0;yi<ci;yi++)
                    {
                        if((l1==Point_2(chull[yi][0],chull[yi][1])))
                            if((l2==Point_2(chull[(yi+1)%ci][0],chull[(yi+1)%ci][1])))
                                fl=0;
                    }
                    if(fl==0)
                    {
                        con[cind][0]=hit->face()->dual()->point();
                        con[cind][1]=hit->opposite()->face()->dual()->point();
                        cind++;
                    }
                }
            }
            else
                if(vd.check_status(hit->target()->point()))
                {
                    l1=hit->face()->dual()->point();
                    l2=hit->opposite()->face()->dual()->point();
                    int fl=1;
                    for(int yi=0;yi<ci;yi++)
                    {
                        if((l1==Point_2(chull[yi][0],chull[yi][1])))
                            if(l2==Point_2(chull[(yi+1)%ci][0],chull[(yi+1)%ci][1]))
                                fl=0;
                    }
                    if(fl==0)
                    {
                        con[cind][0]=hit->face()->dual()->point();
                        con[cind][1]=hit->opposite()->face()->dual()->point();
                        cind++;
                    }
                }
        }
        hit++;
    }while(hit!=vd.halfedges_end());
	
	
    VD::Halfedge_iterator hit3=vd.halfedges_begin();
    do{
        if(hit3->has_source()&&hit3->has_target())
        {
            if(vd.check_status(hit3->source()->point())==0&&vd.check_status(hit3->target()->point())==0)
            {
                emat[emi][0]=hit3->source()->point();
                emat[emi][1]=hit3->target()->point();
                emi++;
            }
        }
        hit3++;
    }while(hit3!=vd.halfedges_end());

 hit=vd.halfedges_begin();
    mi=0;
    do{
        if(hit->has_source()&&hit->has_target())
        {
            if((vd.check_status(hit->source()->point())==1&&vd.check_status(hit->target()->point())==1))
            {
                mat[mi][0]=hit->source()->point();
                mat[mi][1]=hit->target()->point();
                mi++;
            }
        }
        hit++;
    }while(hit!=vd.halfedges_end());
    }
complete:
    if(allconvex==1)
    {
        for(int j=0;j<ci;j++)
        {
            con[j][0]=Point_2(chull[j][0],chull[j][1]);
            con[j][1]=Point_2(chull[(j+1)%ci][0],chull[(j+1)%ci][1]);
        }
        cind=ci;
    }
      
   
    
  collectBoundaryIndices(ci, cind, pcon, con);
}

VoronoiCurve::VoronoiCurve(vector<pair<double, double> > pointVec)
{
        points=pointVec;
	reconstruct();
		
}

vector<pair<int,int>> *VoronoiCurve::getBoundary()
{
	return &_boundary;
}
}

