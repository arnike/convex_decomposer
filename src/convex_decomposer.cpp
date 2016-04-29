#include <string.h>
#include <stdio.h>

#include "hacdCircularList.h"
#include "hacdVector.h"
#include "hacdICHull.h"
#include "hacdGraph.h"
#include "hacdHACD.h"

#include <BulletCollision/CollisionShapes/btCompoundShape.h>
#include <LinearMath/btTransform.h>
#include <LinearMath/btVector3.h>
#include <LinearMath/btAlignedObjectArray.h>

#include "cd_types.h"
#include "cd_wavefront.h"

//sEnableSAT creates the data structures required for performing SAT tests between convex polyhedra, as alternative to GJK
bool sEnableSAT = false;

int main(int argc, char** argv)
{
  //
  // 1. Loading the mesh
  //
  ConvexDecomposition::WavefrontObj wo;
  unsigned int tcount = wo.loadObj("input.obj");

  if (tcount == 0)
  {
    printf("WavefrontObj the object is empty\n");
    return 1;
  }

  btVector3 localScaling(6.f, 6.f, 6.f);

  /*
     btTriangleMesh* trimesh = new btTriangleMesh();


     for (int i = 0; i < wo.mTriCount; i++)
     {
     int index0 = wo.mIndices[i*3];
     int index1 = wo.mIndices[i*3+1];
     int index2 = wo.mIndices[i*3+2];

     btVector3 vertex0(wo.mVertices[index0*3], wo.mVertices[index0*3+1],wo.mVertices[index0*3+2]);
     btVector3 vertex1(wo.mVertices[index1*3], wo.mVertices[index1*3+1],wo.mVertices[index1*3+2]);
     btVector3 vertex2(wo.mVertices[index2*3], wo.mVertices[index2*3+1],wo.mVertices[index2*3+2]);

     vertex0 *= localScaling;
     vertex1 *= localScaling;
     vertex2 *= localScaling;

     trimesh->addTriangle(vertex0, vertex1, vertex2);
     }

  // convex variant
  {
  btConvexHullShape* convexShape = new btConvexHullShape();

  // create a hull approximation
  btConvexShape* tmpConvexShape = new btConvexTriangleMeshShape(trimesh);

  btShapeHull* hull = new btShapeHull(tmpConvexShape);
  btScalar margin = tmpConvexShape->getMargin();
  hull->buildHull(margin);
  tmpConvexShape->setUserPointer(hull);

  printf("# triangles: %d -> %d\n", hwo.mTriCount, ull->numTriangles());
  printf("#   indices: %d -> %d\n", wo.mTriCount*3, hull->numIndices());
  printf("#  vertices: %d -> %d\n", wo.mVertexCount, hull->numVertices());
  printf("reducing vertices by creating a convex hull\n");

  bool updateLocalAabb = false;
  for (int i = 0; i < hull->numVertices(); i++)
  convexShape->addPoint(hull->getVertexPointer()[i], updateLocalAabb);	
  convexShape->recalcLocalAabb();

  if (sEnableSAT)
  convexShape->initializePolyhedralFeatures();

  delete tmpConvexShape;
  delete hull;

  m_collisionShapes.push_back(convexShape);

  float mass = 1.f;

  btTransform startTransform;
  startTransform.setIdentity();
  startTransform.setOrigin(btVector3(0,2,14));

  localCreateRigidBody(mass, startTransform, convexShape);
  }

  // fixed, non-moving concave variant
  {
  bool useQuantization = true;
  btCollisionShape* concaveShape = new btBvhTriangleMeshShape(trimesh, useQuantization);
  startTransform.setOrigin(convexDecompositionObjectOffset);
  localCreateRigidBody(0.f, startTransform, concaveShape);
  m_collisionShapes.push_back(concaveShape);
  }*/

  //-----------------------------------
  // Bullet Convex Decomposition
  //-----------------------------------

  FILE* outputFile = fopen("output.obj", "wb");

  const unsigned int depth = 5;
  const float cpercent     = 5;
  const float ppercent     = 15;
  const unsigned int maxv  = 16;
  const float skinWidth    = 0.0;

  printf("WavefrontObj num triangles read %i\n", tcount);

  ConvexDecomposition::DecompDesc desc;
  desc.mVcount       = wo.mVertexCount;
  desc.mVertices     = wo.mVertices;
  desc.mTcount       = wo.mTriCount;
  desc.mIndices      = (unsigned int *)wo.mIndices;
  desc.mDepth        = depth;
  desc.mCpercent     = cpercent;
  desc.mPpercent     = ppercent;
  desc.mMaxVertices  = maxv;
  desc.mSkinWidth    = skinWidth;

  //-----------------------------------------------
  // HACD
  //-----------------------------------------------

  std::vector<HACD::Vec3<HACD::Real> > points;
  std::vector<HACD::Vec3<long> > triangles;

  for(int i = 0; i < wo.mVertexCount; i++) 
  {
    int index = i*3;
    HACD::Vec3<HACD::Real> vertex(wo.mVertices[index], wo.mVertices[index+1], wo.mVertices[index+2]);
    points.push_back(vertex);
  }

  for(int i = 0; i < wo.mTriCount; i++)
  {
    int index = i*3;
    HACD::Vec3<long> triangle(wo.mIndices[index], wo.mIndices[index+1], wo.mIndices[index+2]);
    triangles.push_back(triangle);
  }

  HACD::HACD hacd;
  hacd.SetPoints(&points[0]);
  hacd.SetNPoints(points.size());
  hacd.SetTriangles(&triangles[0]);
  hacd.SetNTriangles(triangles.size());
  hacd.SetCompacityWeight(0.1);
  hacd.SetVolumeWeight(0.0);

  // HACD parameters
  // Recommended parameters: 2 100 0 0 0 0
  size_t nClusters = 2;
  double concavity = 100;
  bool invert = false;
  bool addExtraDistPoints = false;
  bool addNeighboursDistPoints = false;
  bool addFacesPoints = false;       

  hacd.SetNClusters(nClusters);                     // minimum number of clusters
  hacd.SetNVerticesPerCH(100);                      // max of 100 vertices per convex-hull
  hacd.SetConcavity(concavity);                     // maximum concavity
  hacd.SetAddExtraDistPoints(addExtraDistPoints);   
  hacd.SetAddNeighboursDistPoints(addNeighboursDistPoints);   
  hacd.SetAddFacesPoints(addFacesPoints); 

  hacd.Compute();
  nClusters = hacd.GetNClusters();  

  hacd.Save("output.wrl", false);

  btCompoundShape* compound = new btCompoundShape();
  //m_collisionShapes.push_back(compound);

  btTransform trans;
  trans.setIdentity();

  btVector3	centroid(0, 0, 0);

  size_t base_count = 0;
  for (int c = 0; c < nClusters; c++)
  {
    //generate convex result
    size_t nPoints = hacd.GetNPointsCH(c);
    size_t nTriangles = hacd.GetNTrianglesCH(c);

    float* vertices = new float[nPoints*3];
    unsigned int* triangles = new unsigned int[nTriangles*3];

    HACD::Vec3<HACD::Real> * pointsCH = new HACD::Vec3<HACD::Real>[nPoints];
    HACD::Vec3<long> * trianglesCH = new HACD::Vec3<long>[nTriangles];
    hacd.GetCH(c, pointsCH, trianglesCH);

    // points
    for(size_t v = 0; v < nPoints; v++)
    {
      vertices[3*v] = pointsCH[v].X();
      vertices[3*v+1] = pointsCH[v].Y();
      vertices[3*v+2] = pointsCH[v].Z();
    }

    // triangles
    for(size_t f = 0; f < nTriangles; f++)
    {
      triangles[3*f] = trianglesCH[f].X();
      triangles[3*f+1] = trianglesCH[f].Y();
      triangles[3*f+2] = trianglesCH[f].Z();
    }

    delete [] pointsCH;
    delete [] trianglesCH;

    // Saving the result
    {
      fprintf(outputFile, "## Hull with %d vertices and %d triangles.\r\n", nPoints, nTriangles);
      fprintf(outputFile, "usemtl Material%i\r\n", base_count);
      fprintf(outputFile, "o Object%i\r\n", base_count);

      //calc centroid, to shift vertices around center of mass
      centroid.setValue(0, 0, 0);
      {
        //const unsigned int *src = result.mHullIndices;
        for (unsigned int i = 0; i < nPoints; i++)
        {
          const float *p = &vertices[i*3];
          fprintf(outputFile, "v %0.9f %0.9f %0.9f\r\n", p[0], p[1], p[2]);
          btVector3 vertex(p[0], p[1], p[2]);
          vertex *= localScaling;
          centroid += vertex;
        }
        centroid *= 1.f/(float(nPoints));
      }

      {
        const unsigned int *src = triangles;
        for (unsigned int i = 0; i < nTriangles; i++)
        {
          unsigned int index0 = *src++;
          unsigned int index1 = *src++;
          unsigned int index2 = *src++;

          btVector3 vertex0(vertices[index0*3], vertices[index0*3+1], vertices[index0*3+2]);
          btVector3 vertex1(vertices[index1*3], vertices[index1*3+1], vertices[index1*3+2]);
          btVector3 vertex2(vertices[index2*3], vertices[index2*3+1], vertices[index2*3+2]);

          vertex0 *= localScaling;
          vertex1 *= localScaling;
          vertex2 *= localScaling;

          vertex0 -= centroid;
          vertex1 -= centroid;
          vertex2 -= centroid;

          //trimesh->addTriangle(vertex0, vertex1, vertex2);

          index0 += base_count;
          index1 += base_count;
          index2 += base_count;

          fprintf(outputFile, "f %d %d %d\r\n", index0 + 1, index1 + 1, index2 + 1);
        }
      }
    }
  }

  fclose(outputFile);
}
