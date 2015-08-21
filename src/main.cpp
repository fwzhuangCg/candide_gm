//
//  main.cpp
//  candide
//
//  Created by Damian Stewart on 03.11.13.
//  Copyright (c) 2013 Damian Stewart. All rights reserved.
//

#include <iostream>
#include <unistd.h>

#include "CDApp.h"

#include "CDAssimpLoader.h"
#include "CDMeanValueMeshDeformer.h"
#include <glm/glm.hpp>

using namespace glm;

//************************ASM***************************
#include <vector>
#include <string>

#include <cv.h>
#include <highgui.h>

#include "asmfitting.h"
#include "vjfacedetect.h"
#include "video_camera.h"

#include "levmarmirror.h"

using namespace std;

static void usage_fit()
{
	printf("Usage: fit -m model_file -h cascade_file "
		"{-i image_file | -v video_file | -c camera_idx} -n n_iteration\n\n\n");
	exit(0);
}


static void DrawResult(IplImage* image, const asm_shape& shape)
{
	for(int j = 0; j < shape.NPoints(); j++)
		cvCircle(image, cvPoint(shape[j].x, shape[j].y), 2, CV_RGB(255, 0, 0));
}
//************************ASM***************************


int runCandideApp( int argc, const char* argv[] )
{
	try {
		
		CDApp candideApp(argc, argv);
		
		return candideApp.run();
		
	} catch ( const std::exception &e ) {
		
		std::cerr << "Unhandled exception: "<< e.what() << std::endl;
	}
	
	return 0;
}




int main(int argc, const char * argv[])
{
//*****************ASM**************************
	asmfitting fit_asm;
	const char* model_name = NULL;
	const char* cascade_name = NULL;
	const char* filename = NULL;
	int use_camera = 0;
	int image_or_video = -1;
	int i;
	int n_iteration = 20;
	int camera_idx = 0;

	if(1 == argc)	usage_fit();
	for(i = 1; i < argc; i++)
	{
		if(argv[i][0] != '-') usage_fit();
		if(++i > argc)	usage_fit();
		switch(argv[i-1][1])
		{
		case 'm':
			model_name = argv[i];
			break;
		case 'h':
			cascade_name = argv[i];
			break;
		case 'i':
			if(image_or_video >= 0 || use_camera)
			{
				fprintf(stderr, "only process image/video/camera once\n");
				usage_fit();
			}
			filename = argv[i];
			image_or_video = 'i';
			break;
		case 'v':
			if(image_or_video >= 0 || use_camera)
			{
				fprintf(stderr, "only process image/video/camera once\n");
				usage_fit();
			}
			filename = argv[i];
			image_or_video = 'v';
			break;
		case 'c':
			if(image_or_video >= 0)
			{
				fprintf(stderr, "only process image/video/camera once\n");
				usage_fit();
			}
			use_camera = 1;
			camera_idx = atoi(argv[i]);
			break;
		case 'H':
			usage_fit();
			break;
		case 'n':
			n_iteration = atoi(argv[i]);
			break;
		default:
			fprintf(stderr, "unknown options\n");
			usage_fit();
		}
	}

	double t = (double)cvGetTickCount();
	if(fit_asm.Read(model_name) == false)
		return -1;
	t = ((double)cvGetTickCount() -  t )/  (cvGetTickFrequency()*1000.);
	printf("ASM model file read time cost: %.2f millisec\n", t);
	
	t = (double)cvGetTickCount();
	if(init_detect_cascade(cascade_name) == false)
		return -1;
	t = ((double)cvGetTickCount() -  t )/  (cvGetTickFrequency()*1000.);
	printf("Opencv haar-like file read time cost: %.2f millisec\n", t);
	
	// case 1: process video, here we assume that the video contains only one face,
	// if not, we process with the most central face
	if(image_or_video == 'v')
	{
		int frame_count;
		asm_shape shape, detshape;
		bool flag = false;
		IplImage* image; 
		/* NOTE: the image must not be released, it will be dellocated automatically
		by the class asm_cam_or_avi*/
		int j;
	
		frame_count = open_video(filename);
		if(frame_count == -1)	return false;

		cvNamedWindow("ASM-Search",1);

		for(j = 0; j < frame_count; j ++)
		{
			double t = (double)cvGetTickCount();
			printf("Tracking frame %04i: ", j);
			
			image = read_from_video(j);
			if(image == NULL) continue;
			
			if(j == 0 || flag == false)
			{
				//Firstly, we detect face by using Viola_jones haarlike-detector
				flag = detect_one_face(detshape, image);
				
				//Secondly, we initialize shape from the detected box
				if(flag) 
				{	
					InitShapeFromDetBox(shape, detshape, fit_asm.GetMappingDetShape(), fit_asm.GetMeanFaceWidth());
				}
				else goto show;
			}
			
			//Thirdly, we do image alignment 
			flag = fit_asm.ASMSeqSearch(shape, image, j, true, n_iteration);
			for( int i=0; i<68; ++i ){       //gzc 
				std::cout<<shape[i].x<<" , "<<shape[i].y<<std::endl;    //gzc
			}        //gzc

			//If success, we draw and show its result
			cvWaitKey();
			if(flag) DrawResult(image, shape);
show:
			cvShowImage("ASM-Search", image);
			cvWaitKey(1);
			
			t = ((double)cvGetTickCount() -  t )/  (cvGetTickFrequency()*1000.);
			printf("ASM fitting time cost: %.2f millisec\n", t);
		}

		close_video();
	}
	// case 2: process image, we can process multi-person image alignment
	// also you can process single face alignment by coding like this
	// asm_shape detshape, shape;	
	// bool flag = face_detect.DetectCentralFace(detshape, image);
	// if(flag) asm_common::InitShapeFromDetBox(shape, detshape, 
	//		fit_asm.GetMappingDetShape(), fit_asm.GetMeanFaceWidth());
	// fit_asm.Fitting(shape, image, n_iteration);
	// shape.Write(stdout); //print result
	// for(int l = 0; l < shape.NPoints(); l++)
	//		printf("(%g, %g) ", shape[i].x, shape[i].y);
	else if(image_or_video == 'i')
	{
		IplImage * image = cvLoadImage(filename, 1);
		if(image == 0)
		{
			fprintf(stderr, "Can not Open image %s\n", filename);
			exit(0);
		}

		double t = (double)cvGetTickCount();
		int nFaces;
		asm_shape *shapes = NULL, *detshapes = NULL;
		
		// step 1: detect face
		bool flag =detect_all_faces(&detshapes, nFaces, image);
		t = ((double)cvGetTickCount() -  t )/  (cvGetTickFrequency()*1000.);
		printf("Opencv face detect time cost: %.2f millisec\n", t);
		
		// step 2: initialize shape from detect box
		if(flag)
		{
			shapes = new asm_shape[nFaces];
			for(int i = 0; i < nFaces; i++)
			{
				InitShapeFromDetBox(shapes[i], detshapes[i], fit_asm.GetMappingDetShape(), fit_asm.GetMeanFaceWidth());
			}
		}
		else 
		{
			fprintf(stderr, "This image doesnot contain any faces!\n");
			exit(0);
		}
		
		// step 3: image alignment fitting
		t = (double)cvGetTickCount();
		fit_asm.Fitting2(shapes, nFaces, image, n_iteration);
		t = ((double)cvGetTickCount() -  t )/  (cvGetTickFrequency()*1000.);
		printf("ASM fitting time cost: %.2f millisec\n", t);
					
		// step 4: draw and show result in GUI
		for(int i = 0; i < nFaces; i++)
		{
			DrawResult(image, shapes[i]);
		}
		
		cvSaveImage("result.jpg", image);
		cvNamedWindow("Fitting", 1);
		cvShowImage("Fitting", image);	
		cvWaitKey(0);			
		cvReleaseImage(&image);


		int nmarks = shapes->NPoints();
		int marks[200];
		for(int i = 0; i < nmarks; i++)
		{
			marks[2 * i] = (*shapes)[i].x;
			marks[2 * i + 1] = (*shapes)[i].y;
		}

		LevMarMirror LM;
		LM.LevMarMirroring( marks );


		// step 5: free resource
		delete[] shapes;
		free_shape_memeory(&detshapes);
	}
	// case 3: process camera
	else if(use_camera)
	{
		asm_shape shape, detshape;
		bool flag = false;
		IplImage* image; 
		int j = 0;
		
		if(open_camera(camera_idx) == false)
		{
			fprintf(stderr, "Can not open camera [%d]\n", camera_idx);
			exit(0);
		}

		cvNamedWindow("ASM-Search",1);

		while(1)
		{
			// NOTE: when the parameter is set 1, we can read from camera
			image = read_from_camera();
			if(image == NULL) continue;
			
			if(flag == false)
			{
				//Firstly, we detect face by using Viola_jones haarlike-detector
				flag = detect_one_face(detshape, image);
				
				//Secondly, we initialize shape from the detected box
				if(flag)
				{
					InitShapeFromDetBox(shape, detshape, fit_asm.GetMappingDetShape(), fit_asm.GetMeanFaceWidth());
					j ++;
				}
				else 
					goto show2;
			}
			
			//Thirdly, we do image alignment 
			flag = fit_asm.ASMSeqSearch(shape, image, j, true, n_iteration);
			
			//If success, we draw and show its result
			if(flag) DrawResult(image, shape);
show2:
			cvShowImage("ASM-Search", image);
			cvWaitKey(1);
		}

		close_camera();

	}
//********************ASM**********************************


//*****************************CANDIDE****************************
	char buf[4096];
	getcwd(buf, 4096);                                   //gzc : Get the name of the current working directory
	std::cout<<"cwd: "<<buf<<std::endl;
	
	//runMeanValueMeshDeformerTest();
	return runCandideApp( argc, argv ); 

//*****************************CANDIDE****************************
}
