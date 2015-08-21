#include "levmarmirror.h"
#include "CDFaceData.h"

#define TWO_ROUND
#ifndef TWO_ROUND
int LevMarMirror::m1 = 17;
int LevMarMirror::m2 = 0;
#else
int LevMarMirror::m1 = 9;
int LevMarMirror::m2 = 8;
#endif
int LevMarMirror::n = 74;
int LevMarMirror::AsmPoints[] = {14, 0, 11, 3, 10, 4, 8, 6, 7, //9 points
                                 15, 16, 18, 21, 22, 24, 20, 26, 67, 41, 43, 39, //12points
                                 44, 38, 42, 40, 47, 46, 54, 48, 51, 57, 61, 64,
                                 63,65, 62, 60}; //16 points
int LevMarMirror::CandideVertices[] = {29, 62, 28, 61, 30, 63, 32, 65, 10,
                                       15, 16, 17, 48, 49, 50, 18, 51, 5, 6, 26, 59,
                                       92, 93, 111, 112, 75, 76, 31, 64, 7, 8, 40, 87,
                                       81, 82, 83, 84};
double MatchWeights[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                         1.35, 1.35, 1.35, 1.35, 1.35, 1.35, 1.35, 1.35, 1.35,
                         1.35, 1.35,1.35,1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.25, 1.25,
                         1.25, 1.25, 1.25, 1.25,1.1, 1.1, 1.1, 1.1};
#ifndef TWO_ROUND
int LevMarMirror::CandideParams1[] = {0, 1, 2, 3, 4, 6, 7, 12, 23, 9, 13, 14, 15, 17, 18,19, 29};
int LevMarMirror::CandideParams2[];
#else
int LevMarMirror::CandideParams1[] = {0, 1, 2, 3, 4, 6, 7, 12, 23};
int LevMarMirror::CandideParams2[] = {9, 13, 14, 15, 17, 18, 19, 29};
#endif

#ifndef TWO_ROUND
double LevMarMirror::p1[17];
double LevMarMirror::p2[];
#else
double LevMarMirror::p1[9];
double LevMarMirror::p2[8];
#endif

LevMarMirror::LevMarMirror()
{
}

int LevMarMirror::LevMarMirroring(int *features)
{
    double *x = new double[n];
    for(int i = 0; i < m1; i++)
    p1[i] = 0.0;//translate scalar rotation defermaton parm
    for(int i = 0; i < m2; i++)
    p2[i] = 0.0;
    for(int i = 0; i < n; i++)
        x[i] = features[2 * AsmPoints[i / 2] + i % 2] * MatchWeights[i / 2];
    CDFaceData DataModel( "../data/candide3.wfm");
    CDLog<<"finished load Model";
    CDMesh tempModel = DataModel.getCDMesh();
    CDLog<<"finished get mesh ";
    tempModel.setshapeUnits( DataModel.getShapeUnits() );
    CDLog<<"finished  getShapeUnits";
    tempModel.setanimationUnits( DataModel.getAnimationUnits() );
    CDLog<<"finished getAnimationUnits";
    int it1 = dlevmar_dif(candide1, p1, x, m1, n, 500, NULL, NULL, NULL, NULL, (void *)&tempModel);
    int it2 = dlevmar_dif(candide2, p2, x, m2, n, 500, NULL, NULL, NULL, NULL, (void *)&tempModel);
    tempModel.print();
    CDLog<<it1 << it2;
    delete x;
    std::cout<<it1<<" "<<it2<<std::endl;
    return it1 + it2;
}

void LevMarMirror::candide1(double *p, double *x, int m, int n, void *data)
{
    CDLog<<"before LM";
    register int i;
    // eruFace::Model* CandideModel = (eruFace::Model*) data;
    CDMesh* CandideModel = ( CDMesh* )data;
    std::vector<double> v;
    // CandideModel->getAllParams(v);
    CDLog<<"before getAllParams";
    CandideModel->getAllParams(v);
    for(i = 0; i < m; i++)
    {
        v[CandideParams1[i]] = p[i];
        CDLog << p[i];
    }
    CDLog<<"before setAllParams";
    CandideModel->setAllParams(v);
    CDLog<<"before updateImageCoords";
    CandideModel->updateImageCoords(3, 480, 360, ortographic);
    CDLog<<"after updateImageCoords";
    //CandideModel->setAllParams(v);
    //CandideModel->updateImageCoords(3, 480, 360);
    for(i = 0; i < n; i++)
    {
        x[i] = CandideModel->imageCoord(CandideVertices[i / 2])[i % 2] * MatchWeights[i/ 2];
    }
}

void LevMarMirror::candide2(double *p, double *x, int m, int n, void *data)
{
    register int i;
    // eruFace::Model* CandideModel = (eruFace::Model*) data;
    CDMesh* CandideModel = ( CDMesh* )data;
    std::vector<double> v;
    CandideModel->getAllParams(v);
    for(i = 0; i < m; i++)
    {
        v[CandideParams2[i]] = p[i];
    }
    CandideModel->setAllParams(v);
     CandideModel->updateImageCoords(3, 480, 360, ortographic);
    for(i = 0; i < n; i++)
    {
        x[i] = CandideModel->imageCoord(CandideVertices[i / 2])[i % 2] * MatchWeights[i/ 2];
    }
}

