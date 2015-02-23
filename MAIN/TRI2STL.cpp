#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include "geometry.h"
#include "TRI2STL.h"
#include "smooth.h"
#include "Util.h"
#include "trimesh.h"
#include "Spacing_Field.h"
#include "PList.h"
#include "Octree_Storage.h"
#include "tetgen.h"

//global variables
extern FILE *in_f, *jou_f, *out_f, *debug_f; /*input output journal files global*/

void tri2stl_lib(geometry *geom, int nregions, int *regions)
{
  int bdim = 132;
  char *sname = new char[bdim];
  const int bdim2 = 400;
  char buff[bdim2];
  int i, j, k, l, m, n;
  int n0, n1, n2;
  Vector v0, v1, norm;
  FILE *tet_f = NULL;
  FILE *ucd_f = NULL;
  //declared here for use w or w/out restart
  POLYMESH *finalmesh = new POLYMESH();
  POLYMESH *tetmesh = new POLYMESH();
  
  finalmesh->nelems = 0;
  finalmesh->nn = geom->n_gnodes;
  finalmesh->nb = geom->ngb;
  tetmesh->nelems = 0;
  tetmesh->nn = 0;
  tetmesh->nb = 0;
  
  fprintf(out_f,"\nAllocating for geom boundary.\n");
  fflush(out_f);
  
  finalmesh->boundary = new BOUNDARY_MESH[geom->ngb];
  
  //set all names and elists
  for (i = 0; i < geom->ngb; i++)
    {
    finalmesh->boundary[i].name = new char[400];
    sprintf(finalmesh->boundary[i].name,"%s",geom->g_bname[i+1]);
    finalmesh->boundary[i].elist = new List();
    }

  fprintf(out_f,"\nAllocating for geom boundary elements.\n");
  fflush(out_f);
  
  if ((finalmesh->nelems+geom->n_gfacets) > finalmesh->element_dim)
    my_mem(finalmesh->element_dim, finalmesh->nelems, &(finalmesh->element), (finalmesh->nelems+geom->n_gfacets), ELEMENT_CHUNK);
  
  //since all elements are one face, set up nf and f_n now
  for (k = finalmesh->nelems; k < (finalmesh->nelems+geom->n_gfacets); k++)
    {
    finalmesh->element[k].Initialize();
    finalmesh->element[k].nf = 1;
    finalmesh->element[k].f_n = new List*[1];
    finalmesh->element[k].f_n[0] = new List();
    }
    
  fprintf(out_f,"\nStoring facets as boundary objects.\n");
  fflush(out_f);
  
  //loop through geom facets, noting all unique nodes at this point
  //set up geom as elements, then tris as elements
  for (i=1; i <= geom->ngb; i++)
    {
    for (j=geom->n_begin[i]; j < geom->n_end[i]; j++)
      {
      n0 = geom->g_facet[j].nodes[0];
      n1 = geom->g_facet[j].nodes[1];
      n2 = geom->g_facet[j].nodes[2];
      
      //take the resulting triangles and set up tri
      finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(n0);
      finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(n1);
      finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(n2);
      
      finalmesh->boundary[i-1].elist->Add_To_List(finalmesh->nelems); 
      //inc nelems
      finalmesh->nelems++;
      }
    }
  
  //add in geom nodes
  for (i = 0; i < geom->n_gnodes; i++)
    {
    if ((finalmesh->nn + 1) > finalmesh->node_dim)
      my_mem(finalmesh->node_dim, finalmesh->nn, &(finalmesh->node), (finalmesh->nn + 1), NODE_CHUNK);
      
    finalmesh->node[i].vert = geom->g_vert[i];
    }
      
  //finally, output stl file
  sname[0] = '\0';
  
  sprintf(sname,"FASTAR.stl");

  if ((ucd_f=fopen(sname,"w")) == NULL)
    {
    fprintf(stderr,"\nCouldn't open surface mesh output file %s",sname);
    fflush(stderr);
    exit(0);
    }
   
  //STL format
  fprintf(ucd_f, "solid for_pointwise\n");
  fflush(ucd_f);
    
  Vector areasum = Vector(0.0,0.0,0.0);
    
  //put in triangles and normals
  for (i = 0; i < finalmesh->nb; i++)
    {
    for (j = 0; j < finalmesh->boundary[i].elist->max; j++)
      {
      n0 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[0];
      n1 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[1];
      n2 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[2];
      
      //find normal vector
      v0 = Vector(finalmesh->node[n0].vert,finalmesh->node[n1].vert);
      v1 = Vector(finalmesh->node[n0].vert,finalmesh->node[n2].vert);
      norm = v0 % v1;
      norm.normalize();
        
      fprintf(ucd_f, "facet normal %16.10e %16.10e %16.10e\n",norm[0],norm[1],norm[2]);
        
      //do area summation for inner bd only..reset norm in process
      norm = (v0%v1)*0.5;
      areasum += norm;
        
      fprintf(ucd_f, "  outer loop\n");
      fprintf(ucd_f, "    vertex %16.10e %16.10e %16.10e\n",finalmesh->node[n0].vert[0],finalmesh->node[n0].vert[1],finalmesh->node[n0].vert[2]);
      fprintf(ucd_f, "    vertex %16.10e %16.10e %16.10e\n",finalmesh->node[n1].vert[0],finalmesh->node[n1].vert[1],finalmesh->node[n1].vert[2]);
      fprintf(ucd_f, "    vertex %16.10e %16.10e %16.10e\n",finalmesh->node[n2].vert[0],finalmesh->node[n2].vert[1],finalmesh->node[n2].vert[2]);
      fprintf(ucd_f, "  endloop\n");
      fprintf(ucd_f, "endfacet\n");
      }
    }
      
  fprintf(ucd_f, "endsolid for_pointwise\n");
  fflush(ucd_f);
    
  fprintf(out_f, "Area summation (x, y, z) = %16.10e %16.10e %16.10e\n",areasum[0],areasum[1],areasum[2]);
  fflush(out_f);
    
  //close ucd file
  fclose(ucd_f);

  tetgenio *in, *out;
  in = new tetgenio();
  in->initialize();
  out = new tetgenio();
  out->initialize();
  tetgenio::facet *f;
  tetgenio::polygon *p;
  tetgenbehavior flags;

  in->firstnumber = 1;
  in->numberofpoints = finalmesh->nn;

  fprintf(out_f,"\npoints=%d\n",in->numberofpoints);
  fflush(out_f);

  in->pointlist = new REAL[in->numberofpoints * 3];
  for (i = 0; i < finalmesh->nn; i++)
    {
    j = i*3;
    in->pointlist[j] = (REAL) finalmesh->node[i].vert[0];
    in->pointlist[j+1] = (REAL) finalmesh->node[i].vert[1];
    in->pointlist[j+2] = (REAL) finalmesh->node[i].vert[2];
    }

  in->numberoffacets = finalmesh->nelems;

  fprintf(out_f,"\nfacets=%d\n",in->numberoffacets);
  fflush(out_f);    

  in->facetlist = new tetgenio::facet[in->numberoffacets];
  in->facetmarkerlist = new int[in->numberoffacets];
  k = 0; //for numbering
  for (i = 0; i < geom->ngb; i++)
    {
    for (j = 0; j < finalmesh->boundary[i].elist->max; j++)
      {
      n0 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[0];
      n1 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[1];
      n2 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[2];
      
      f = &(in->facetlist[k]);
      in->init(f);
      f->numberofpolygons = 1;
      f->polygonlist = new tetgenio::polygon[1];
      f->numberofholes = 0;
      f->holelist = NULL;
      p = &(f->polygonlist[0]);
      in->init(p);
      p->numberofvertices = 3;
      p->vertexlist = new int[p->numberofvertices];
      p->vertexlist[0] = n0+1;
      p->vertexlist[1] = n1+1;
      p->vertexlist[2] = n2+1;
      in->facetmarkerlist[k] = i;
      k++;
      }
    }
  in->save_nodes("TRI2STL");
  in->save_poly("TRI2STL");  
    
  char *options = new char[100];
  strcpy(options,"pqQAA"); //use extra Y for full recovery, one less for don't care at all!
  flags.parse_commandline(options);

  tetrahedralize(&flags, in, out);

  out->save_nodes("TRI2STLout");
  out->save_elements("TRI2STLout");
  out->save_faces("TRI2STLout");

  //in->deinitialize();
  //delete in;
  //delete [] options; 
  
  //vars for tetmesh
  int tetnn = 0, ntet = 0;
  int ntribd = 0;
  int **tribd = 0;
  int **tets = 0;
  Point *tetnodes = 0;
  double node0, node1, node2;

  //now, read in files output by tet mesher
  if ((tet_f=fopen("TRI2STLout.node","r")) == NULL)
    {
    fprintf(stderr,"\nCouldn't open tet meshfile TRI2STLout.node.  Exiting....");
    fflush(stderr);
    exit(0);
    }
      
  //read num nodes
  fgets(buff,bdim2,tet_f);
  sscanf(buff,"%d %d %d %d",&tetnn,&i,&i,&i);
    
  //now, allocate
  /*fprintf(out_f,"\nAllocating for %d nodes.\n", out->numberofpoints);
  fflush(out_f);
  tetnodes = (Point*)malloc(out->numberofpoints*sizeof(Point));
  tetmesh->node = new NODE[out->numberofpoints];*/
  fprintf(out_f,"\nAllocating for %d nodes.\n", tetnn);
  fflush(out_f);
  tetnodes = (Point*)malloc(tetnn*sizeof(Point));
    
  //create nodes
  for (i = 0; i < tetnn; i++)
    {
    fgets(buff,bdim2,tet_f);
    sscanf(buff,"%d %lg %lg %lg",&j,&(node0),&(node1),&(node2));
    tetnodes[i] = Point(node0,node1,node2);
    }
  /*for (i = 0; i < out->numberofpoints; i++)
    {
    tetnodes[i] = Point(out->pointlist[3*i],out->pointlist[3*i+1],out->pointlist[3*i+2]);
    }

  tetnn = out->numberofpoints;*/

  fclose(tet_f);
  tet_f = NULL;

  //open elem file just to get number of tets
  if ((tet_f=fopen("TRI2STLout.ele","r")) == NULL)
    {
    fprintf(stderr,"\nCouldn't open tet meshfile TRI2STLout.ele.  Exiting....");
    fflush(stderr);
    exit(0);
    }
      
  //read num tets
  fgets(buff,bdim2,tet_f);
  sscanf(buff,"%d %d %d",&ntet,&i,&i);
    
  fprintf(out_f,"\nNumber of tets = %d.\n", ntet);
  fflush(out_f);

  fclose(tet_f);
  tet_f = NULL;
      
  if ((tet_f=fopen("TRI2STLout.face","r")) == NULL)
    {
    fprintf(stderr,"\nCouldn't open tet meshfile TRI2STLout.face.  Exiting....");
    fflush(stderr);
    exit(0);
    }
      
  //read num faces
  fgets(buff,bdim2,tet_f);
  sscanf(buff,"%d %d",&ntribd,&i);
    
  //fprintf(out_f,"\nAllocating for %d facets.\n", out->numberoffacets);
  fprintf(out_f,"\nAllocating for %d facets.\n", ntribd);
  fflush(out_f);

  //tribd = (int**)calloc(out->numberoffacets,sizeof(int*));
  //for (k = 0; k < out->numberoffacets; k++)
  tribd = (int**)calloc(ntribd,sizeof(int*));
  for (k = 0; k < ntribd; k++)
    {
    tribd[k] = (int*)calloc(3,sizeof(int));
    for (l = 0; l < 3; l++)
      tribd[k][l] = -1;
    }

  //ntribd = out->numberoffacets;

  //tetmesh->element = new POLY_ELEMENT[out->numberoffacets + out->numberoftetrahedra];
  tetmesh->element = new POLY_ELEMENT[ntribd + ntet];
  //reset nelems
  tetmesh->nelems = 0;

  //set up bd mesh
  //set number of boundaries...make last one internal bd, go back and delete later
  /*tetmesh->nb = geom->ngb; 
  
  //need to add boundaries
  tetmesh->boundary = new BOUNDARY_MESH[geom->ngb];
    
  //first, set all names
  for (i = 0; i < geom->ngb; i++)
    {
    tetmesh->boundary[i].name = new char[400];
    sprintf(finalmesh->boundary[i].name,"%s",geom->g_bname[i+1]);
    }
     
  //now, set up lists!
  for (i = 0; i < geom->ngb+1; i++)
    {
    tetmesh->boundary[i].elist = new List();
    }*/
    
  tetmesh->nb = geom->ngb;
  tetmesh->boundary = new BOUNDARY_MESH[geom->ngb];
  
  List **bound = new List*[geom->ngb];
  for (i = 0; i < geom->ngb; i++)
    bound[i] = new List();
  
  //set all names and elists
  for (i = 0; i < geom->ngb; i++)
    {
    tetmesh->boundary[i].name = new char[400];
    sprintf(tetmesh->boundary[i].name,"%s",geom->g_bname[i+1]);
    tetmesh->boundary[i].elist = new List();
    }
 
  //for (j = 0; j < out->numberoffacets; j++)
  for (j = 0; j < ntribd; j++)
    {
    /*f = &(out->facetlist[j]);
    p = &(f->polygonlist[0]);
    tribd[j][0] = p->vertexlist[0] - 1;
    tribd[j][1] = p->vertexlist[1] - 1;
    tribd[j][2] = p->vertexlist[2] - 1;*/

    //read faces
    fgets(buff,bdim2,tet_f);
    sscanf(buff,"%d %d %d %d %d",&i,&(tribd[j][0]),&(tribd[j][1]),&(tribd[j][2]),&k);
    
    tribd[j][0]--;
    tribd[j][1]--;
    tribd[j][2]--;
    bound[k]->Add_To_List(j);
    }

  fclose(tet_f);
  tet_f = NULL;

  //open elem file
  if ((tet_f=fopen("TRI2STLout.ele","r")) == NULL)
    {
    fprintf(stderr,"\nCouldn't open tet meshfile TRI2STLout.ele.  Exiting....");
    fflush(stderr);
    exit(0);
    }
      
  //read num tets
  fgets(buff,bdim2,tet_f);
  sscanf(buff,"%d %d %d",&ntet,&i,&i);
    
  //fprintf(out_f,"\nAllocating for %d tets.\n", out->numberoftetrahedra);
  fprintf(out_f,"\nAllocating for %d tets.\n", ntet);
  fflush(out_f);
  //tets = (int**)calloc(out->numberoftetrahedra,sizeof(int*));
  //for (j = 0; j < out->numberoftetrahedra; j++)
  tets = (int**)calloc(ntet,sizeof(int*));
  for (j = 0; j < ntet; j++)
    {
    tets[j] = (int*)calloc(4,sizeof(int));
    for (k = 0; k < 4; k++)
      tets[j][k] = -1;
    } 
    
  //begin reading tets...need only those in region
  int t1, t2, t3, t4, reg, newntet;

  //keeps track of used tets
  newntet = 0;

  //make list of regions
  List reglist;
  reglist.construct();

  for (j = 0; j < nregions; j++)
    reglist.Add_To_List(regions[j]);

  //for (j = 0; j < out->numberoftetrahedra; j++)
  for (j = 0; j < ntet; j++)
    {
    /*tets[j][0] = out->tetrahedronlist[4*j] - 1;
    tets[j][1] = out->tetrahedronlist[4*j+1] - 1;
    tets[j][2] = out->tetrahedronlist[4*j+2] - 1;
    tets[j][3] = out->tetrahedronlist[4*j+3] - 1;*/

    //read faces
    fgets(buff,bdim2,tet_f);
    sscanf(buff,"%d %d %d %d %d %d",&i,&(t1),&(t2),&(t3),&(t4),&(reg));

    if (!reglist.Is_In_List(reg))
        continue;

    tets[newntet][0] = t1-1;
    tets[newntet][1] = t2-1;
    tets[newntet][2] = t3-1;
    tets[newntet][3] = t4-1;

    //inc nelems
    newntet++;
    }
 
  //NOTE: ntet is not nelems...includes tets we dropped and nelems includes bd faces
  reglist.destruct();
  //finally, reset ntet
  ntet = newntet;
    
  fclose(tet_f);
  tet_f = NULL;

  //finally, dump extra nodes and realloc tets
  tets = (int**)realloc((void*)tets,ntet*sizeof(int*));

  //set up map array
  int *map = new int[tetnn];

  //init
  for (i = 0; i < tetnn; i++)
    map[i] = -1;

  for (i = 0; i < ntet; i++)
    for (j = 0; j < 4; j++)
      map[tets[i][j]] = 0;

  for (i = 0; i < ntribd; i++)
    for (j = 0; j < 3; j++)
      map[tribd[i][j]] = 0;

  j = 0;
  for (i = 0; i < tetnn; i++)
    if (map[i] == 0)
      map[i] = j++;

  //remap tets and bd faces
  for (i = 0; i < ntet; i++)
    for (j = 0; j < 4; j++)
      tets[i][j] = map[tets[i][j]];
 
  for (i = 0; i < ntribd; i++)
    for (j = 0; j < 3; j++)
      tribd[i][j] = map[tribd[i][j]];

  j = 0;
  //finally, reset nodes
  for (i = 0; i < tetnn; i++)
    if (map[i] >= 0)
      {
      tetnodes[map[i]] = Point(tetnodes[i][0],tetnodes[i][1],tetnodes[i][2]);
      j++;
      }

  //reset tetnn
  tetnn = j;

  //realloc
  tetnodes = (Point*)realloc((void*)tetnodes,tetnn*sizeof(Point));

  //now, reset ntet
  //ntet = out->numberoftetrahedra;
  
  tetmesh->node = new NODE[tetnn];
  tetmesh->nn = tetnn;
    
  //create nodes
  for (i = 0; i < tetnn; i++)
    {
    tetmesh->node[i].vert = tetnodes[i];
    }
  /*for (i = 0; i < out->numberofpoints; i++)
    {
    tetmesh->node[i].vert = Point(out->pointlist[3*i],out->pointlist[3*i+1],out->pointlist[3*i+2]);
    }

  tetmesh->nn = out->numberofpoints;*/
  
  for (i = 0; i < geom->ngb; i++)
    {
    for (j = 0; j < bound[i]->max; j++)
      {
      //set up as elements
      tetmesh->element[tetmesh->nelems].nf = 1;
      tetmesh->element[tetmesh->nelems].f_n = new List*[1];
      tetmesh->element[tetmesh->nelems].f_n[0] = new List();

      tetmesh->element[tetmesh->nelems].f_n[0]->Add_To_List(tribd[bound[i]->list[j]][0]);
      tetmesh->element[tetmesh->nelems].f_n[0]->Add_To_List(tribd[bound[i]->list[j]][1]);
      tetmesh->element[tetmesh->nelems].f_n[0]->Add_To_List(tribd[bound[i]->list[j]][2]);

      //tetmesh->boundary[out->facetmarkerlist[j]].elist->Add_To_List(tetmesh->nelems);
      tetmesh->boundary[i].elist->Add_To_List(tetmesh->nelems);
      
      //inc nelems
      tetmesh->nelems++;
      }
    }
    
  //for (j = 0; j < out->numberoftetrahedra; j++)
  for (j = 0; j < ntet; j++)
    {
    //set up as elements
    tetmesh->element[tetmesh->nelems].nf = 4;
    tetmesh->element[tetmesh->nelems].f_n = new List*[4];
    for (k = 0; k < 4; k++)
      tetmesh->element[tetmesh->nelems].f_n[k] = new List();

    tetmesh->element[tetmesh->nelems].f_n[0]->Add_To_List(tets[j][0]);
    tetmesh->element[tetmesh->nelems].f_n[0]->Add_To_List(tets[j][2]);
    tetmesh->element[tetmesh->nelems].f_n[0]->Add_To_List(tets[j][1]);

    tetmesh->element[tetmesh->nelems].f_n[1]->Add_To_List(tets[j][0]);
    tetmesh->element[tetmesh->nelems].f_n[1]->Add_To_List(tets[j][1]);
    tetmesh->element[tetmesh->nelems].f_n[1]->Add_To_List(tets[j][3]);

    tetmesh->element[tetmesh->nelems].f_n[2]->Add_To_List(tets[j][1]);
    tetmesh->element[tetmesh->nelems].f_n[2]->Add_To_List(tets[j][2]);
    tetmesh->element[tetmesh->nelems].f_n[2]->Add_To_List(tets[j][3]);

    tetmesh->element[tetmesh->nelems].f_n[3]->Add_To_List(tets[j][2]);
    tetmesh->element[tetmesh->nelems].f_n[3]->Add_To_List(tets[j][0]);
    tetmesh->element[tetmesh->nelems].f_n[3]->Add_To_List(tets[j][3]);
      
    //inc nelems
    tetmesh->nelems++;
    }

  //out->deinitialize();
  //delete out;

  //finally write out tetmesh and delete finalmesh
  //reset file name before volume mesh output
  sname[0] = '\0';
  
  sprintf(sname,"TRI2STL_tet.cgns");

  //output vol mesh
  tetmesh->smooth_io(1,geom,sname);
  
  fprintf(out_f,"\nFinished with output tet CGNS file.\n");
  fflush(out_f);

  delete [] map;
  for (i = 0; i < geom->ngb; i++)
    delete bound[i];
  delete bound;

  delete finalmesh;
  delete tetmesh;
  
    
  //finally, delete stitching mem
  free(tetnodes);
  for (i = 0; i < ntet; i++)
    free(tets[i]);
  free(tets);
  for (j = 0; j < ntribd; j++)
    {
    free(tribd[j]);
    }
  free(tribd);
      
  return;
}
