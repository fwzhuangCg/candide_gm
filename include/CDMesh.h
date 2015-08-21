//
//  CDMesh.h
//  candide
//
//  Created by damian on 03/11/13.
//  Copyright (c) 2013 Damian Stewart. All rights reserved.
//

#ifndef __candide__CDMesh__
#define __candide__CDMesh__

#include <iostream>
#include <vector>
#include <glm/glm.hpp>
#include <Fl/gl.h>
#include <map>
//**********testmodel_0***************
#include <unordered_map>
//**********testmodel_0***************
//***********************************
//class CDFaceDistortionUnit;
class CDMesh;
enum Projection { ortographic, weakPerspective, perspective };
class CDMeshOperation
{
public:
	/*! @abstract Transform vertices in targetMesh by transform; rebuild normals if necessary. */
	static CDMesh transform( const CDMesh& mesh, const glm::mat4 transform );
};


class CDFaceDistortionUnit
{
public:
	CDFaceDistortionUnit( const std::string& name = "" );
	void setName( const std::string& n ) { name = n; }
	std::string getName() const { return name; }
	
	void addVertexOffset( unsigned short vIdx, glm::vec3 offset );
	
	/// apply this shape unit to the target mesh with the given factor (usually 0..1)
	void apply( float amount, CDMesh& target ) const;

	std::map<unsigned short, glm::vec3> offsets;
	
private:
	
	std::string name;
	
};

class CDMesh
{
	friend class CDFaceDistortionUnit;
	
public:
	CDMesh()
	{
		 _rotation.x = 0.0;
		 _rotation.y = 0.0;
		 _rotation.z = 0.0;
		 _scale.x=  1.0;
		 _scale.y = 1.0; 
		 _scale.z=  1.0; 
		 _translation.x = 0.0;
		 _translation.y = 0.0;
		 _translation.z = 0.0;
		 _staticParams.assign(nStaticDeformations(), 0.0);
		 _dynamicParams.assign(nDynamicDeformations(), 0.0);
	};
	virtual ~CDMesh() {};
	
	void clear() { vertices.clear(); triangles.clear(); }
	
	void addVertex( const glm::vec3 &v ) { vertices.push_back(v); }
	void addVertexNormal( const glm::vec3 &v ) { vertexNormals.push_back(v); }
	void addTextureCoordinate( const glm::vec2 &coord ) { textureCoordinates.push_back(coord); }
	void addFace( int v0, int v1, int v2 );
	
	void draw() const;
	void drawBoundingBox() const;
	
	struct Triangle {
		GLushort v[3];
		Triangle(): Triangle(0,0,0) {};
		Triangle( short v0, short v1, short v2 ) { v[0] = v0; v[1] = v1; v[2] = v2; }
	};
	
	size_t getNumVertices() const { return vertices.size(); }
	const glm::vec3& getVertex(size_t which) const { return vertices.at(which); }
	void setVertex(size_t which, const glm::vec3& v) { vertices[which] = v; }
	
	size_t getNumTriangles() const { return triangles.size(); }
	const Triangle& getTriangle(size_t which) const { return triangles.at(which); }
	
	size_t getNumNormals() const { return vertexNormals.size(); }
	const glm::vec3& getNormal(size_t which) const { return vertexNormals.at(which); }
	
	size_t getNumTextureCoordinates() const { return textureCoordinates.size(); }
	const glm::vec2& getTextureCoordinate( size_t i ) const { return textureCoordinates.at(i); }
	void removeAllTextureCoordinates() { textureCoordinates.clear(); }
	
	glm::vec3 getBoundingBoxCenter() const;
	glm::vec3 getBoundingBoxSize() const;
	void getBoundingBox( glm::vec3& minCornerOut, glm::vec3& maxCornerOut ) const;

	void updateNormals();
	//add zfw 

	inline glm::vec3& imageCoord( int v ) { return _imageCoords[v]; }

            	//**************************testmodel_1*******************************

	int   nStaticDeformations() const;
    	int   nDynamicDeformations() const;

   	void  setStaticParam( int, double );
    	double getStaticParam( int ) const;

   	 void  setDynamicParam( int, double );
    	double  getDynamicParam( int ) const;

	void   setAllParams( const std::vector<double>& params );
	void   getAllParams( std::vector<double>& );	
	//**************************testmodel_1*******************************

	void  vertexscale ( std::vector<glm::vec3> & vs , const glm::vec3 & v) ;
	void inline vertexscale (std::vector<glm::vec3> & vs,  double f )  { glm::vec3 v(f,f,f); vertexscale(vs , v); }
	void inline vertexscale( std::vector<glm::vec3> & vs, double x, double y, double z = 1 ) { glm::vec3 v(x,y,z); vertexscale(vs, v); }

	void  vertextranslate(std::vector<glm::vec3> & vs, const  glm::vec3 &v );
	void inline vertextranslate( std::vector<glm::vec3> & vs, double x, double y, double z = 0 ) { glm::vec3  v(x,y,z); vertextranslate(vs, v); }
	void  vertexrotate(std::vector<glm::vec3> & vs, const  glm::vec3& );
	void  vertexrotate(std::vector<glm::vec3> & vs, double, double, double ) ;

	void   vertextransform( std::vector<glm::vec3> & vs, const glm::mat3& );

	void   vertexapplyDeformation(std::vector<glm::vec3>&, const CDFaceDistortionUnit& , double );
	void   applyDeformations  ( std::vector<glm::vec3>&, const std::vector<CDFaceDistortionUnit>&, const std::vector<double>& );
	bool  updateStatic();
            	bool  updateDynamic();
	void  updateGlobal();
	void   updateImageCoords( int viewPortSize, int imageWidth, int imageHeight, Projection projection );

	void  setshapeUnits( const std::vector<CDFaceDistortionUnit> &SU ){ shapeUnits = SU; }
	void  setanimationUnits( const std::vector<CDFaceDistortionUnit> &AU ){ animationUnits = AU; }

	void print() const
	{
		std::cout<< _transformedCoords.size()<<std::endl;
		for(auto itr = _transformedCoords.begin(); itr != _transformedCoords.end(); itr++)
		{
			std::cout<<(*itr)[0] << "    "<<(*itr)[1] << "    "<<(*itr)[2] << std::endl;
		}
	}


	//add zfw
	
protected:

	virtual void setupArrays() const;
	virtual void teardownArrays() const;
	//*************************testmodel_2********************************
	std::vector<double>  _staticParams;
	std::unordered_map<std::string, int> _staticIndices;
	bool _staticParamsModified;

	std::vector<double>  _dynamicParams;
	std::unordered_map<std::string, int> _dynamicIndices;
	bool _dynamicParamsModified;

	//*************************testmodel_2********************************
	
private:
	
	
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> vertexNormals;
	std::vector<glm::vec2> textureCoordinates;

	std::vector<Triangle> triangles;

	//add zfw
	
	glm::vec3  _rotation;
	glm::vec3  _scale; 
	glm::vec3  _translation;
	glm::mat3   _transform;
	


	std::vector<CDFaceDistortionUnit> shapeUnits;
	std::vector<CDFaceDistortionUnit> animationUnits;

	//add zfw
	std::vector<glm::vec3>  _imageCoords;
	std::vector<glm::vec3> _staticCoords;
	std::vector<glm::vec3> _dynamicCoords;
	std::vector<glm::vec3> _transformedCoords;

};

#endif /* defined(__candide__CDMesh__) */

