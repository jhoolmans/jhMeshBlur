//
//  File: jhMeshBlur.cpp
//
//  Author: Jeroen Hoolmans (jhoolmans@gmail.com)
//
//  Please share any changes you make to this code and mention my name if you use or share this code.
//

#include <iostream>

#include <maya/MPxDeformerNode.h> 

#include <maya/MFnMesh.h>
#include <maya/MMatrix.h>

#include <maya/MFnPlugin.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MFnData.h>

#include <maya/MTypeId.h> 
#include <maya/MPlug.h>

#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MDagPath.h>
#include <maya/MArrayDataHandle.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnUnitAttribute.h>

#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MGlobal.h>
#include <maya/MTime.h>
#include <maya/MAnimControl.h>
#include <maya/MItGeometry.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnPointArrayData.h>

#include <maya/MDagModifier.h>
#include <maya/MFloatVectorArray.h>
#include <maya/MVectorArray.h>
#include <maya/MVector.h>

class jhMeshBlur : public MPxDeformerNode
{
public:
                            jhMeshBlur();
    virtual                 ~jhMeshBlur();

    static      void*		creator();
    static      MStatus		initialize();

    // deformation function
    virtual     MStatus     deform(MDataBlock& block,  MItGeometry& iter, const MMatrix& mat, unsigned int multiIndex);

    // local node attributes
    static      MObject     aOldMeshData;       // save old mesh inside this attribute
    static      MObject		aMeshInput;         // the basic input
    static      MObject     aStrength;          // strength of the blur
    static      MObject     aTreshhold;         // make small distance changes not affecting the blur
    static      MObject     aShapeFactor;       // TODO: Make this have influence over the shape, changing the position of the center mass
    static      MObject     aTime;              // Needed when doing the Cache etc
    static      MObject     aTweakBlur;         // turn tweakmode on and off
    static      MObject     aQuadInterp;        // Switch interpolation mode between Linear and Quad
    static      MObject     aInterpPower;       // Strength of the interpolation
    static      MTypeId		id;                 // the node id (default)

private:
};


MTypeId		 jhMeshBlur::id( 0x00117840 );

// local attributes
//
MObject     jhMeshBlur::aMeshInput;
MObject     jhMeshBlur::aOldMeshData;
MObject     jhMeshBlur::aStrength;
MObject     jhMeshBlur::aTreshhold;
MObject     jhMeshBlur::aShapeFactor;
MObject     jhMeshBlur::aTime;
MObject     jhMeshBlur::aTweakBlur;
MObject     jhMeshBlur::aQuadInterp;
MObject     jhMeshBlur::aInterpPower;


jhMeshBlur::jhMeshBlur() {}
jhMeshBlur::~jhMeshBlur() {}

void* jhMeshBlur::creator()
{
	return new jhMeshBlur();
}

//
//  Initialize the node
//
MStatus jhMeshBlur::initialize()
{
    // attribute types
    MFnUnitAttribute    unitAttr;
    MFnNumericAttribute	nAttr;
    MFnTypedAttribute   tAttr;

    aOldMeshData = tAttr.create("oldMesh","om",MFnData::kPointArray);
    tAttr.setArray(true);
    tAttr.setHidden(true);
    tAttr.setIndexMatters(true);

    // create the attributes
    aStrength = nAttr.create( "Strength", "str", MFnNumericData::kFloat,1.0);
    nAttr.setStorable(true);
    nAttr.setKeyable(true);
    nAttr.setMax(1.0);
    nAttr.setMin(0.0);

    aTreshhold = nAttr.create( "Treshold", "tres", MFnNumericData::kFloat,0.0);
    nAttr.setStorable(true);
    nAttr.setKeyable(true);
    nAttr.setMin(0.0);

    aShapeFactor = nAttr.create( "ShapeFactor", "shapef", MFnNumericData::kFloat,0.5);
    nAttr.setStorable(true);
    nAttr.setKeyable(true);
    nAttr.setMax(1.0);
    nAttr.setMin(0.0);

    aTweakBlur = nAttr.create( "TweakBlur", "tweak", MFnNumericData::kBoolean,false);
    nAttr.setKeyable(false);
    nAttr.setChannelBox(true);

    aQuadInterp = nAttr.create( "QuadInterpolation", "qi", MFnNumericData::kBoolean,true);
    nAttr.setKeyable(false);
    nAttr.setChannelBox(true);

    aInterpPower = nAttr.create( "InterpolationPower", "interp", MFnNumericData::kDouble, 0.75);
    nAttr.setKeyable(true);
    nAttr.setMax(1.0);
    nAttr.setMin(0.0);

    aTime = unitAttr.create( "time", "tm", MFnUnitAttribute::kTime, 1.0 );
    unitAttr.setStorable(true);
    unitAttr.setCached(true);
	unitAttr.setReadable(true);
	unitAttr.setWritable(true);
	unitAttr.setAffectsAppearance(true);
	unitAttr.setAffectsWorldSpace(true);

    // Make the attributes visible to the user
    addAttribute( aStrength);
    addAttribute( aTreshhold);
    addAttribute( aTime);
    addAttribute( aTweakBlur);
    addAttribute( aQuadInterp);
    addAttribute( aInterpPower);
    addAttribute( aOldMeshData);

    // Make sure when an attribute changes, the node updates
    attributeAffects( aTime, outputGeom );
	attributeAffects( aStrength, outputGeom );
    attributeAffects( aTreshhold, outputGeom );
    attributeAffects( aQuadInterp, outputGeom );
    attributeAffects( aInterpPower, outputGeom );

    // Not implented yet, but make the weights paintable :)
    MGlobal::executeCommand("makePaintable -attrType multiFloat -sm deformer jhMeshBlur weights;");

    return MStatus::kSuccess;
}

//
//      Linear interpolation between 0 to 1
//
double linearInterp(double given, double a, double b){
    return (given - a) * (1.0 / (b - a));
}

//
//      Quadratic curves According to Wikipedia
//
MPoint quadInterpBetween(MPoint a, MPoint b, MPoint c, double t){
    MPoint EdgeA = a + ((b-a) * t);
    MPoint EdgeB = b + ((c-b) * t);
	
    return EdgeA + (EdgeB - EdgeA) * t;
}

//
//      Deform computation
//
MStatus jhMeshBlur::deform( MDataBlock& block,MItGeometry& iter,const MMatrix& m,unsigned int multiIndex)
{
    MStatus returnStatus;

    // Envelope
    float envData = block.inputValue(envelope, &returnStatus).asFloat();
	CHECK_MSTATUS(returnStatus);

    if(envData == 0)
		return MS::kFailure;

    /*
     VARIABLES
     */
    //float factor = block.inputValue(aShapeFactor, &returnStatus).asFloat();
    float fStrength = block.inputValue(aStrength, &returnStatus).asFloat();
	CHECK_MSTATUS(returnStatus);
	
	if (fStrength == 0)
		return MS::kFailure;
	
    float fThreshold = block.inputValue(aTreshhold, &returnStatus).asFloat();
	CHECK_MSTATUS(returnStatus);
    float fW = 0.0f; // weight
    float fDistance;
    fStrength *= envData;

    double dKracht = block.inputValue(aInterpPower, &returnStatus).asDouble();
	CHECK_MSTATUS(returnStatus);
    double dDotProduct;  // Dotproduct of the point

    bool bTweakblur = block.inputValue(aTweakBlur, &returnStatus).asBool();
	CHECK_MSTATUS(returnStatus);
	
    bool bQuad = block.inputValue(aQuadInterp, &returnStatus).asBool();
	CHECK_MSTATUS(returnStatus);
	
	MTime inTime = block.inputValue(aTime).asTime();
    int nTijd = (int)inTime.as(MTime::kFilm);


    MFloatVectorArray currentNormals;   // normals of mesh
    MFnPointArrayData fnPoints;         // help converting to MPointArrays
    MFloatVector dirVector;             // direction vector of the point
    MFloatVector normal;                // normal of the point
    MPointArray savedPoints;            // save all point before edited
    MMatrix matInv = m.inverse();       // inversed matrix
    MPoint ptA;                         // current point (iter mesh)
    MPoint ptB;                         // previous point (iter mesh)
    MPoint ptC;                         // mesh before previous point (iter mesh)

    // get node, use node to get inputGeom, use inputGeom to get mesh data, use mesh data to get normal data
    MFnDependencyNode nodeFn(this->thisMObject());

    MPlug inGeomPlug(nodeFn.findPlug(this->inputGeom,true));
    MObject inputObject(inGeomPlug.asMObject());
    MFnMesh inMesh(inputObject);

    inMesh.getVertexNormals(true, currentNormals);

    // get the previous mesh data
    MPlug oldMeshPlug = nodeFn.findPlug(MString("oldMesh"));
    MPlug oldMeshPositionsAPlug = oldMeshPlug.elementByLogicalIndex((multiIndex*4) + 0);
    MPlug oldMeshPositionsBPlug = oldMeshPlug.elementByLogicalIndex((multiIndex*4) + 1);
    MPlug oldMeshPositionsCPlug = oldMeshPlug.elementByLogicalIndex((multiIndex*4) + 2); // cache for tweak mode
    MPlug oldMeshPositionsDPlug = oldMeshPlug.elementByLogicalIndex((multiIndex*4) + 3); // cache for tweak mode

    // convert to MPointArrays
    MObject objOldMeshA;
    MObject objOldMeshB;
    MObject objOldMeshC; // cache
    MObject objOldMeshD; // cache

    oldMeshPositionsAPlug.getValue(objOldMeshA);
    oldMeshPositionsBPlug.getValue(objOldMeshB);
    oldMeshPositionsCPlug.getValue(objOldMeshC); // cache
    oldMeshPositionsDPlug.getValue(objOldMeshD); // cache

    fnPoints.setObject(objOldMeshA);
    MPointArray oldMeshPositionsA = fnPoints.array();
    
    fnPoints.setObject(objOldMeshB);
    MPointArray oldMeshPositionsB = fnPoints.array();
    
    fnPoints.setObject(objOldMeshC);
    MPointArray oldMeshPositionsC = fnPoints.array(); // cache
    
    fnPoints.setObject(objOldMeshD);
    MPointArray oldMeshPositionsD = fnPoints.array(); // cache

    
    
    // If mesh position variables are empty,fill them with default values
    if(oldMeshPositionsA.length() == 0 || nTijd <= 1){
        iter.allPositions(oldMeshPositionsA);

        for(int i=0; i < oldMeshPositionsA.length(); i++)
        {
            // convert to world
            oldMeshPositionsA[i] = oldMeshPositionsA[i] * m;
        }
		
        oldMeshPositionsB.copy(oldMeshPositionsA);
        oldMeshPositionsC.copy(oldMeshPositionsA); // cache
        oldMeshPositionsD.copy(oldMeshPositionsA); // cache
    }
	
	// get back old date again
	if (bTweakblur == true) { // restore cache
		oldMeshPositionsA.copy(oldMeshPositionsC);
		oldMeshPositionsB.copy(oldMeshPositionsD);
	}
    
    
    iter.allPositions(savedPoints);
    for(int i=0; i < savedPoints.length(); i++)
    {
        // convert points to world points
        savedPoints[i] = savedPoints[i] * m;
    }

    // Actual Iteration through points
    for (; !iter.isDone(); iter.next()){
        // get current position
        ptA = iter.position();
        // get old positions
        ptB = oldMeshPositionsA[iter.index()] * matInv;
        ptC = oldMeshPositionsB[iter.index()] * matInv;

        fDistance = ptA.distanceTo(ptB);
        fW = weightValue(block,multiIndex,iter.index());


        if (fDistance * (fStrength*fW) < fThreshold && fThreshold > 0){
            iter.setPosition(ptA);
        } else {
            // aim/direction vector to calculate strength
            dirVector = (ptA - ptB); // (per punt)
            dirVector.normalize();

            normal = currentNormals[iter.index()];

            dDotProduct = normal.x * dirVector.x + normal.y * dirVector.y + normal.z * dirVector.z;

            
            if(bQuad == true){
                MVector vecA(((ptB - ptC) + (ptA - ptB)) / 2);
                vecA.normalize();

                MPoint hiddenPt(ptB + (vecA * fDistance) * dKracht);
                ptA = quadInterpBetween(ptB, hiddenPt, ptA, (1 - fStrength * fW) + (linearInterp(dDotProduct, -1, 1) * (fStrength * fW) ) );
            } else {
                MPoint halfway = (ptA - ptB) * 0.5;
                MPoint offset = halfway * dDotProduct * (fStrength*fW);
                ptA = ptA - ((halfway * (fStrength*fW)) - offset); // + (offset * strength);
            }
            // set new value

            iter.setPosition(ptA);
        }
    }
    if(bTweakblur == false){
        oldMeshPositionsD.copy(oldMeshPositionsB);
        oldMeshPositionsC.copy(oldMeshPositionsA);
        oldMeshPositionsB.copy(oldMeshPositionsA);
        oldMeshPositionsA.copy(savedPoints);

        // Save back to plugs
        objOldMeshA = fnPoints.create(oldMeshPositionsA);
        objOldMeshB = fnPoints.create(oldMeshPositionsB);
        objOldMeshC = fnPoints.create(oldMeshPositionsC);
        objOldMeshD = fnPoints.create(oldMeshPositionsD);
		
        oldMeshPositionsAPlug.setValue(objOldMeshA);
        oldMeshPositionsBPlug.setValue(objOldMeshB);
        oldMeshPositionsCPlug.setValue(objOldMeshC);
        oldMeshPositionsDPlug.setValue(objOldMeshD);
    }
    
    return returnStatus;
}





// standard initialization procedures
//

MStatus initializePlugin( MObject obj )
{
	MStatus result;
	MFnPlugin plugin( obj, "Jeroen Hoolmans", "0.66", "Any");
	result = plugin.registerNode( "jhMeshBlur", jhMeshBlur::id, jhMeshBlur::creator, jhMeshBlur::initialize, MPxNode::kDeformerNode );

	return result;
}

MStatus uninitializePlugin( MObject obj)
{
	MStatus result;
	MFnPlugin plugin( obj );
	result = plugin.deregisterNode( jhMeshBlur::id );
	return result;
}
