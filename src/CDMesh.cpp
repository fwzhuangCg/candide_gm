//
//  CDMesh.cpp
//  candide
//
//  Created by damian on 03/11/13.
//  Copyright (c) 2013 Damian Stewart. All rights reserved.
//

#include "CDMesh.h"
#include "CDUtilities.h"

#include <Fl/gl.h>

using namespace glm;
using namespace std;


void CDMesh::addFace( int v0, int v1, int v2 )
{
	assert(v0<vertices.size() && "v0 out of bounds");
	assert(v1<vertices.size() && "v1 out of bounds");
	assert(v2<vertices.size() && "v2 out of bounds");
	
	triangles.push_back(Triangle(v0,v1,v2));
	
}


void CDMesh::setupArrays() const
{
	
	glEnableClientState( GL_VERTEX_ARRAY );
	glVertexPointer(3, GL_FLOAT, 0, &vertices[0][0] );
	
	bool doNormals = (vertexNormals.size()==vertices.size());
	if ( doNormals ) {
		glEnableClientState(GL_NORMAL_ARRAY);
		glEnable(GL_NORMALIZE);
		glNormalPointer(GL_FLOAT, 0, &vertexNormals[0][0] );
		glEnable( GL_LIGHTING );
	} else {
		glDisable( GL_LIGHTING );
	}
	
	if ( textureCoordinates.size() == getNumVertices() ) {
		glEnableClientState( GL_TEXTURE_COORD_ARRAY );
		glTexCoordPointer( 2, GL_FLOAT, 0, &textureCoordinates[0][0] );
	}
}


void CDMesh::draw() const
{
	setupArrays();

	glDrawElements( GL_TRIANGLES, (GLsizei)triangles.size()*3, GL_UNSIGNED_SHORT, &triangles[0].v[0]);
	
	teardownArrays();
	
}

void CDMesh::teardownArrays() const
{
	glDisableClientState( GL_VERTEX_ARRAY );
	glDisableClientState( GL_TEXTURE_COORD_ARRAY );
	
	bool doNormals = (vertexNormals.size()==vertices.size());
	if ( doNormals ) {
		glDisableClientState(GL_NORMAL_ARRAY);
		glDisable( GL_LIGHTING );
	}
}

void CDMesh::drawBoundingBox() const
{
	vec3 vertices[8];
	
	vec3 centre = getBoundingBoxCenter();
	vec3 bounds = getBoundingBoxSize();
	vec3 halfBounds = bounds*0.5f;
	for ( int i=0;i<2;i++ )
		for ( int j=0;j<2;j++ )
			for ( int k=0;k<2;k++ )
				vertices[i*4+j*2+k] = vec3(centre.x+halfBounds.x*(i==0?-1:1),centre.y+halfBounds.y*(j==0?-1:1),centre.z+halfBounds.z*(k==0?-1:1));
	
	short lines[12*2] = { 0,1, 0,4, 4,5, 5,1, 2,3, 2,6, 6,7, 3,7, 0,2, 1,3, 4,6, 5,7 };
	glEnableClientState( GL_VERTEX_ARRAY );
	glVertexPointer(3, GL_FLOAT, 0, vertices);
	glDrawElements( GL_LINES, 12*2, GL_UNSIGNED_SHORT, &lines[0] );
	
	glDisableClientState(GL_VERTEX_ARRAY);
	
}


void CDMesh::getBoundingBox( vec3& minCorner, vec3& maxCorner ) const
{
	for ( int i=0; i<vertices.size(); i++ ) {
		const vec3& v = vertices[i];
		if ( i==0 ) {
			minCorner = v;
			maxCorner = v;
		} else {
			minCorner.x = std::min(v.x,minCorner.x);
			minCorner.y = std::min(v.y,minCorner.y);
			minCorner.z = std::min(v.z,minCorner.z);
			maxCorner.x = std::max(v.x,maxCorner.x);
			maxCorner.y = std::max(v.y,maxCorner.y);
			maxCorner.z = std::max(v.z,maxCorner.z);
		}
	}
}

vec3 CDMesh::getBoundingBoxCenter() const
{
	vec3 minCorner, maxCorner;
	getBoundingBox(minCorner,maxCorner);
	return (minCorner+maxCorner)*0.5f;
}

vec3 CDMesh::getBoundingBoxSize() const
{
	vec3 minCorner, maxCorner;
	getBoundingBox(minCorner,maxCorner);
	return maxCorner-minCorner;
}

void CDMesh::updateNormals()
{
	vector<pair<vec3,int> > accumulatedVertexNormals;
	accumulatedVertexNormals.resize(vertices.size());
	
	// accumulate vertex normals
	for ( auto t: triangles ) {
		vec3 p0 = vertices[t.v[0]];
		vec3 p1 = vertices[t.v[1]];
		vec3 p2 = vertices[t.v[2]];
		
		vec3 normal = normalize(cross(p1-p0, p1-p2));

		for ( int i=0; i<3; i++ ) {
			auto& vn = accumulatedVertexNormals.at(t.v[i]);
			vn.first += normal;
			vn.second++;
		}
	}
	
	// apply
	vertexNormals.resize(vertices.size());
	for ( int i=0; i<accumulatedVertexNormals.size(); i++ )
	{
		if ( accumulatedVertexNormals[i].second>0 ) {
			vertexNormals.at(i) = normalize(accumulatedVertexNormals[i].first/(float)accumulatedVertexNormals[i].second);
		} else {
			vertexNormals.at(i) = vec3(0,0,1);
		}
	}
	
	
	
}

/// add zfw


//*****************************testmodel_0*************************************
int CDMesh::nStaticDeformations() const
{
	//return _staticDeformations.size();
	return 38;
}

int CDMesh::nDynamicDeformations() const
{
	//return _dynamicDeformations.size();
	return 70;
}

void  CDMesh::setStaticParam(int paramNo, double paramVal)
{
	if (paramNo >= nStaticDeformations() || paramNo < 0)
	{
	    throw std::runtime_error("Static parameter number out of range!");
	}
	_staticParams[paramNo] = paramVal;
	_staticParamsModified = true;
}

double CDMesh::getStaticParam(int paramNo) const
{
	if (paramNo >= nStaticDeformations() || paramNo < 0)
	{
	    throw std::runtime_error("Static parameter number out of range!");
	}
	return _staticParams[paramNo];
}

void CDMesh::setDynamicParam(int paramNo, double paramVal)
{
	if (paramNo >= nDynamicDeformations() || paramNo < 0)
	{
	    throw std::runtime_error("Dynamic parameter number out of range!");
	}
	_dynamicParams[paramNo] = paramVal;
	_dynamicParamsModified = true;
}

double CDMesh::getDynamicParam(int paramNo) const
{
	if (paramNo >= nDynamicDeformations() || paramNo < 0)
	{
	    throw std::runtime_error("Dynamic parameter number out of range!");
	}
	return _dynamicParams[paramNo];
}

void  CDMesh::setAllParams( const std::vector<double>& params )
{
	if (params.size() != 9 + nStaticDeformations() + nDynamicDeformations())
	    {
		throw std::runtime_error("Bad argument: Wrong input vector size!");
	    }
		_rotation[0]    = params[0];
		_rotation[1]    = params[1];
		_rotation[2]    = params[2];
		_scale[0]       = params[3];
		_scale[1]       = params[4];
		_scale[2]       = params[5];
		_translation[0] = params[6];
		_translation[1] = params[7];
		_translation[2] = params[8];
	    for (int p = 0; p < nStaticDeformations(); p++)
	    {
	        setStaticParam(p, params[p+9]);
	    }
	    for (int p = 0; p < nDynamicDeformations(); p++)
	    {
	        setDynamicParam(p, params[p+9+nStaticDeformations()]);
	    }
}

void  CDMesh::getAllParams( std::vector<double>& params)
{
	params.resize(9 + nStaticDeformations() + nDynamicDeformations());
	params[0] = _rotation[0];
	params[1] = _rotation[1];
	params[2] = _rotation[2];
	params[3] = _scale[0];
	params[4] = _scale[1];
	params[5] = _scale[2];
	params[6] = _translation[0];
	params[7] = _translation[1];
	params[8] = _translation[2];

	for (int p = 0; p < nStaticDeformations(); p++)
	{
	    params[p+9] = getStaticParam(p);
	}
	for (int p = 0; p < nDynamicDeformations(); p++)
	{
	    params[p+9+nStaticDeformations()] = getDynamicParam(p);
	}
}


void CDMesh::vertextranslate( std::vector<glm::vec3> & vs, const  glm::vec3 &v )
{
	for(std::vector<glm::vec3>::iterator itr = vs.begin(); itr != vs.end(); itr++)
	{
		*itr += v;	
	}
}

void  CDMesh::vertexscale (std::vector<glm::vec3> & vs,  const glm::vec3 &  v ) 
{
	for(std::vector<glm::vec3>::iterator itr = vs.begin(); itr != vs.end(); itr++)
	{
		*itr  *= v;	
	}
}

void CDMesh::vertexrotate( std::vector<glm::vec3> & vs, const glm::vec3& v )
{
	vertexrotate(vs, v[0], v[1], v[2] );
}

void   CDMesh::vertextransform(std::vector<glm::vec3> & vs,  const glm::mat3& m)
{
	for (int i = 0; i < vs.size(); i++)
	 {
	 	vs[i][0] =  m[0][0]*_transformedCoords[i][0] +  m[0][1]*vs[i][1] +  m[0][2]*vs[i][2];
	 	vs[i][1] =  m[1][0]*_transformedCoords[i][0] +  m[1][1]*vs[i][1] +  m[1][2]*vs[i][2];
	 	vs[i][2] =  m[2][0]*_transformedCoords[i][0] +  m[2][1]*vs[i][1] +  m[2][2]*vs[i][2];
	 }
}

void CDMesh::vertexrotate(std::vector<glm::vec3> & vs, double x, double y, double z )
{
    glm::mat3 m;
    m[0][0] =  cos(z)*cos(y);
    m[0][1] = -sin(z)*cos(x)-cos(z)*sin(y)*sin(x);
    m[0][2] =  sin(z)*sin(x)-cos(z)*sin(y)*cos(x);

    m[1][0] =  sin(z)*cos(y);
    m[1][1] =  cos(z)*cos(x)-sin(z)*sin(y)*sin(x);
    m[1][2] = -cos(z)*sin(x)-sin(z)*sin(y)*cos(x);
 
    m[2][0] =  sin(y);
    m[2][1] =  cos(y)*sin(x);
    m[2][2] =  cos(y)*cos(x);

    vertextransform( vs, m );
}

void CDMesh::vertexapplyDeformation(std::vector<glm::vec3>  & v, const CDFaceDistortionUnit& du, double coeff )
{

	std::map<unsigned short, glm::vec3>::const_iterator  iter;
	for(iter = du.offsets.begin(); iter != du.offsets.end(); iter++)
	{
		int  vertexNo = iter->first;
		v[vertexNo][0]  += coeff * iter->second[0];
		v[vertexNo][1]  += coeff * iter->second[1];
		v[vertexNo][2]  += coeff * iter->second[2];
	}
}

void CDMesh::applyDeformations( std::vector<glm::vec3>  & v, const std::vector<CDFaceDistortionUnit>& dv, const std::vector<double>& p )
{

	for ( std::vector<CDFaceDistortionUnit>::size_type i = 0; i < dv.size(); i++ )  //?changed from i<=0
	{ 
		vertexapplyDeformation(v, dv[i], p[i] );
	}

}

bool CDMesh::updateStatic()
{
	_staticCoords  = vertices;
 	applyDeformations(_staticCoords, shapeUnits, _staticParams);
	_staticParamsModified = false;
	_dynamicParamsModified = true;
	return true;
}

bool CDMesh::updateDynamic()
{
	if (!updateStatic() && !_dynamicParamsModified)
    	{
		return false;
   	 }
   	 _dynamicCoords =  _staticCoords;
    	applyDeformations(_dynamicCoords, animationUnits, _dynamicParams);
	_dynamicParamsModified = false;
	return true;
}

void  CDMesh::updateGlobal()
{
	updateDynamic();
    	_transformedCoords = _dynamicCoords;
	vertexscale(_transformedCoords, _scale);
	vertexrotate(_transformedCoords, _rotation);
	vertextranslate(_transformedCoords, _translation);
}

void  CDMesh::updateImageCoords( int viewPortSize, int imageWidth, int imageHeight, Projection projection)
{
	updateGlobal();
	_imageCoords = _transformedCoords;
	 double f = double(std::max(imageWidth, imageHeight)) / double(viewPortSize);
	switch(projection)
	{
		case ortographic:
	 	{
			vertexscale(_imageCoords, f, -f);
			break;
	 	}
		case weakPerspective:
		 {
		 	 f /= (1.0 + _translation[2]);
			vertexscale(_imageCoords, f, -f);
			break;
		}
		case perspective:
		{
			for (int i = 0; i < _imageCoords.size(); i++)
		            {
				_imageCoords[i] *= (f / (1.0 + _translation[2] + _transformedCoords[i][2]));
			}
			break;
		}
		
	}
	vertextranslate(_imageCoords, double(imageWidth) / 2.0, double(imageHeight) / 2.0);
}


//// add zfw


CDMesh CDMeshOperation::transform( const CDMesh& inputMesh, const glm::mat4 transform )
{
	CDMesh targetMesh = inputMesh;
	
	// transform vertices
	for ( size_t i=0; i<targetMesh.getNumVertices(); i++ ) {
		vec4 v = vec4(targetMesh.getVertex(i),1.0f);
		v = transform*v;
		targetMesh.setVertex(i, vec3(v.x,v.y,v.z) );
	}
	// recalculate normals
	if ( targetMesh.getNumNormals() ) {
		targetMesh.updateNormals();
	}
	
	return targetMesh;
}
