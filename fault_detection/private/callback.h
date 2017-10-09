

/*80*************************************************************************************/
/*80*************************************************************************************/
/*80********************************** NOT TO MODIFY ************************************/
/*80*************************************************************************************/
/*80*************************************************************************************/


#if !defined (LIBRARY_CALLBACK)
	#define LIBRARY_CALLBACK

#if !defined(EXPORT_LIB)
	#define EXPORT_LIB __declspec (dllexport)
#endif

typedef void				(*ndfunc)											();
typedef long				(*ndfunc_function_callback)	  (long id, void * data_event, void * data_callback);
typedef ndfunc			(*attfunc_functionCallback)		(char *);

#define RET_VOID		(void)(( void(*)
#define RET_LONG		(long)(( long(*)
#define RET_PVOID		(void *)(( void *(*)
#define RET_PCHAR		(char *)(( char *(*)
#define ATT_CALL		)(functionCallback

#if defined(LIBRARY_MAIN)	
	attfunc_functionCallback					functionCallback;
#else
	extern attfunc_functionCallback		functionCallback;
#endif



/*80*************************************************************************************/	
//
// Types
//
/*80*************************************************************************************/

// Donnée (type_data)
#define TYPE_BYTE							1				// entier 8 bits non signé [0,255]
#define TYPE_CHAR							2				// entier 8 bits signé [-128,127]
#define TYPE_DOUBLE						3				// virgule flottante sur 64 bits
#define TYPE_FLOAT						4				// virgule flottante sur 32 bits
#define TYPE_SHORT						5				// entier 16 bits signé [-32768,32767]
#define TYPE_LONG							6				// entier 32 bits signé
#define TYPE_RGBA							7				// quadruplet (R,G,B,A) 4x8 bits
#define TYPE_UCHAR						1				// entier 8 bits non signé [0,255]
#define TYPE_USHORT						8				// entier 16 bits non signé [0,65535]
#define TYPE_ULONG						9				// entier 32 bits non signé
	
// Classe (type_class)
#define CLASS_DEFAULT					0
#define CLASS_LABEL						1
#define CLASS_DIRECTION				2
#define CLASS_BOOLEAN					3
#define CLASS_ORIENTATION			4

// Fichier (type_file)
#define DECL_FILE							1
#define DECL_FILE_TRANSPOSE		2	

// Opération de mise à jour sur fichier (type_update)
#define ATT_OP_COPY						0
#define ATT_OP_ADD						1
#define ATT_OP_MUL						2
#define ATT_OP_MIN						3
#define ATT_OP_MAX						4
#define ATT_OP_NONE						99


/*80*************************************************************************************/	
//
// Entrées : déclaration
//
/*80*************************************************************************************/

long ndImageGetIn();
long ndBlockGetIn();
long ndStudyGetIn();
long ndView3dGetIn();


/*80*************************************************************************************/	
//
// Sorties : déclaration, initialisation et création
//
/*80*************************************************************************************/

long ndImageNew();
long ndBlockNew();
long ndGraphNew();
long ndView3dNew();

void ndImageSetTypeData(long id, long type_data);
void ndImageSetSize(long id, long width, long height);
void ndBlockSetTypeData(long id, long type_data);
void ndBlockSetSize(long id, long width, long height, long depth);
void ndBlockSetFile(long id, long type_file, char * filename);
void ndBlockSetMemory(long id);
void ndGraphSetSize(long id, long width, long height);

void ndImageShow(long id);
void ndBlockShow(long id);
void ndGraphShow(long id);
void ndView3dShow(long id);


/*80*************************************************************************************/	
//
// Images 2D
//
/*80*************************************************************************************/

void ndImageGetSize(long id, long * width, long * height);
long ndImageGetClickNumber(long id);
void ndImageGetClickArray(long id, long ** array_x, long ** array_y);
long ndImageGetCaptureNumber(long id);
void ndImageGetCaptureArray(long id, long ** array_x, long ** array_y);
void * ndImageGetData(long id, long type_data);
void * ndImageGetPtrData(long id);
long ndImageGetTypeData(long id);
void ndImageSetClass(long id, long type_class);
void ndImageAddObject2d(long id, long id_object2d);
void ndImageDeleteObject2d(long id, long id_object2d);
void ndImageDeleteObject2dByName(long id, char * name);
long ndImageAddLine(long id, double x1, double y1, double x2, double y2);
long ndImageAddSquare(long id, double x, double y, double c);
long ndImageSaveBitmap(long id, char * filename);


/*80*************************************************************************************/	
//
// Images 3D
//
/*80*************************************************************************************/

void ndBlockGetSize(long id, long * width, long * height, long * depth);
void * ndBlockGetFront(long id, long no, long type_data);
void * ndBlockGetPartialFront(long id, long no, long type_data, long x, long y, long dx, long dy);
void * ndBlockGetRight(long id, long no, long type_data);
void * ndBlockGetPartialRight(long id, long no, long type_data, long x, long y, long dx, long dy);
void * ndBlockGetTop(long id, long no, long type_data);
void * ndBlockGetPartialTop(long id, long no, long type_data, long x, long y, long dx, long dy);
void * ndBlockGetImageFromTopPatch(long id, long no, long type_data, long * patch_dx, long * patch_dy, long nb, long * width, long * height);
void * ndBlockGetImageSONE(long id, long no, long type_data, long * width, long * height);
void * ndBlockGetImageNOSE(long id, long no, long type_data, long * width, long * height);
void ndBlockSetUpdateFile(long id, long type_update);
void ndBlockSetPartialFront(long id, long no, long x, long y, void * data, long w, long h);
void ndBlockSetFront(long id, long no, void * data);
void ndBlockSetFrontByType(long id, long no, void * data, long type_data);
void ndBlockSetPartialFrontByType(long id, long no, long x, long y, void * data, long w, long h, long type_data);
void ndBlockSetPartialRight(long id, long no, long x, long y, void * data, long w, long h);
void ndBlockSetRight(long id, long no, void * data);
void ndBlockSetRightByType(long id, long no, void * data, long type_data);
void ndBlockSetPartialRightByType(long id, long no, long x, long y, void * data, long w, long h, long type_data);
void ndBlockSetPartialTop(long id, long no, long x, long y, void * data, long w, long h);
void ndBlockSetTop(long id, long no, void * data);
void ndBlockSetTopByType(long id, long no, void * data, long type_data);
void ndBlockSetPartialTopByType(long id, long no, long x, long y, void * data, long w, long h, long type_data);
void ndBlockSetImageSONE(long id, long no, void * data);
void ndBlockSetImageSONEByType(long id, long no, void * data, long type_data);
void ndBlockSetImageNOSE(long id, long no, void * data);
void ndBlockSetImageNOSEByType(long id, long no, void * data, long type_data);
long ndBlockGetClickNumber(long id);
void ndBlockGetClickArray(long id, long ** array_x, long ** array_y, long ** array_z);
long ndBlockGetCaptureNumber(long id);
void ndBlockGetCaptureArray(long id, long ** array_x, long ** array_y, long ** array_z);
long ndBlockGetCursor(long id, long * x, long * y, long * z);
void * ndBlockGetPtrData(long id);
long ndBlockGetTypeData(long id);
void ndBlockSetClass(long id, long type_class);
void ndBlockCopyToFile(long id, char * filename);


/*80*************************************************************************************/	
//
// Graphes
//
/*80*************************************************************************************/


/*80*************************************************************************************/	
//
// Vues 3D
//
/*80*************************************************************************************/

void ndView3dGetSize(long id, long * width, long * height, long * depth);
void ndView3dSetSize(long id, long width, long height, long depth);
void ndView3dAddObject3d(long id, long id_object3d);
void ndView3dDeleteObject3d(long id, long id_object3d);
void ndView3dDeleteObject3dByName(long id, char * name);


/*80*************************************************************************************/	
//
// Vidéos
//
/*80*************************************************************************************/

long ndVideoGetId();
void ndVideoSetFormatIn(long id, long format);
void ndVideoSetFormatOut(long id, long format);
void * ndVideoGetDataIn(long id);
void * ndVideoGetDataOut(long id);
long ndVideoGetWidth(long id);
long ndVideoGetHeight(long id);
void * ndVideoGetUserData(long id);
void ndVideoSetUserData(long id, void * user_data);
void ndVideoStart(long id, ndfunc videoOn, ndfunc videoOff);
long ndVideoGetDuration();
long ndVideoGetCurrentDuration();
long ndVideoGetCurrentTime();
long ndVideoGetNumberOfFrames(long id);
long ndVideoGetNoCurrentFrame(long id);
void ndVideoMoveTo(long id, long no);
void ndVideoRefreshOut(long id);


void * ndVideoFileCreate(char * filename);
void ndVideoFileSet(void * avi, long width, long height, double frequency);
void ndVideoFileAddFrame(void * avi, unsigned char * img_rgba, long width, long height);
void * ndVideoFileClose(void * avi);


/*80*************************************************************************************/	
//
// Fenêtres : images 2D, images 3D, vues 3D, vidéos, ...
//
/*80*************************************************************************************/

void ndWindowRedraw(long id_object);
void ndWindowClose(long id_object);
void ndWindowSetName(long id_object, const char * name);
void ndWindowUserDataNew(long id_object, char * name);
void ndWindowUserDataDelete(long id_object, char * name);
void ndWindowUserDataSet(long id_object, char * name, void * ptr);
void * ndWindowUserDataGet(long id_object, char * name);
/* Callbacks (Images 2D)
			"LBUTTONDOWN"				-> "x", "y", "x_image", "y_image", "x_screen", "y_screen"
			"LBUTTONUP"					-> "x", "y", "x_image", "y_image", "x_screen", "y_screen"
			"LBUTTONDBLCLICK"
			"MBUTTONDOWN"				-> "x", "y", "x_image", "y_image", "x_screen", "y_screen"
			"MBUTTONUP"					-> "x", "y", "x_image", "y_image", "x_screen", "y_screen"
			"MBUTTONDBLCLICK"
			"RBUTTONDOWN"				-> "x", "y", "x_image", "y_image", "x_screen", "y_screen"
			"RBUTTONUP"					-> "x", "y", "x_image", "y_image", "x_screen", "y_screen"
			"RBUTTONDBLCLICK"
			"MOUSEMOVE"					-> "x", "y", "x_image", "y_image", "x_screen", "y_screen"	
			"IMAGE_UPDATE"			-> "type_view", "x", "y", "z"
			"MSG_HOOK"					-> "hWnd", "msg", "wparam", "lparam"
*/

void ndWindowCallbackAdd(long id_object, char * name, ndfunc_function_callback function_callback, void * data_callback);
void ndWindowCallbackSub(long id_object, char * name, ndfunc_function_callback function_callback);
long ndWindowCallbackEventGet(long id_object, void * data_event, char * name_event, void * data_to_get);
long ndWindowCallbackEventSet(long id_object, void * data_event, char * name_event, void * data_to_set);
void * ndWindowGetHwnd(long id_object);


/*80*************************************************************************************/	
//
// Divers 
//
/*80*************************************************************************************/

// Etat d'avancement des calculs
void ndProgressSetMinMax(long mini, long maxi);
void ndProgressSetPosition(long position);
void ndProgressReset();

// Fenêtres Popup
long ndPopupImage(void * data, long w, long h, long type_data);
long ndPopupHistogram(void * data, long length, long type_data);

// Libération de ressources
void ndFree(void * ptr);

// Messages
void ndMsg(char * message);

#define ndMsgFormat \
RET_VOID(char * format, ...) \
ATT_CALL("ndMsgFormat")))

// Menus
long ndMenuCreate(char * title);
void ndMenuAddItem(long menu, char * title, ndfunc function);
void ndMenuAddSeparator(long menu);
void ndMenuAddSubMenu(long menu, long submenu);

// Boite de paramètres

#define ndBoxParameter \
RET_LONG(char * string, ...) \
ATT_CALL("ndBoxParameter")))

#define ndBoxParameterAdd \
RET_PVOID(void * data, char * string, ...) \
ATT_CALL("ndBoxParameterAdd")))

long ndBoxParameterDialog(void * data);

// Boite de messages
void ndMsgBox(char * title, char * message);
long ndMsgBoxYNC(char * title, char * message);

// Librairie
void ndLibraryAddMenu(long id_library, char * insert, long menu);
void ndLibraryAddMenuPopup(long id_library, long id_object, long x, long y, long menu);

// Chargement automatique
void ndLoadFromFile(char * filename);

void ndObject3dSetTransparency(long id, double alpha);



/*80*************************************************************************************/	
//
// Objets 2D
//
/*80*************************************************************************************/

typedef struct _PT2D
{
	double x;
	double y;	
} PT2D;
	
#define OBJECT2D_SETOFPOINTS			0
#define OBJECT2D_LINE							1
#define OBJECT2D_POLYGON					2
#define OBJECT2D_CIRCLE						3
#define OBJECT2D_ASTERISK					4
#define OBJECT2D_ARROW						5
#define OBJECT2D_SQUARE						6

long ndObject2dCreate(long type);
void ndObject2dSetName(long id_object2d, char * name);
void ndObject2dSetColor(long id_object2d, double r, double g, double b);
void ndObject2dSetThickness(long id_object2d, long thickness);
void ndObject2dSetArrayPoints(long id_object2d, PT2D * array_pt2d, long nb_pt2d);
void ndObject2dTranslate(long id_object2d, double sx, double sy);
void ndObject2dScale(long id_object2d, double sx, double sy);


/*80*************************************************************************************/	
//
// Objets 3D
//
/*80*************************************************************************************/

typedef struct _PT3D
{
	double x;
	double y;
	double z;
} PT3D;

typedef struct _MESH3D
{
	long nb_pt;
	PT3D * array_pt3d;
	long nb_triangle;
	long * array_triangle;
	void * unknown;
} MESH3D;

void * ndMesh3dCube();
void * ndMesh3dSphere(long precision);
void * ndMesh3dCylinder(long precision);
void * ndMesh3dCone(long precision);
void * ndMesh3dArrow(PT3D * pt_begin, PT3D * pt_end, double radius);
void * ndMesh3dLine(PT3D * array_pt, long nb_pt, double radius, long precision);
void * ndMesh3dReadFromFile(char * filename);
void * ndMesh3dAlloc(long nb_pt, long nb_triangle);
void * ndMesh3dFree(MESH3D * mesh3d);
void ndMesh3dScale(MESH3D * mesh3d, double sx, double sy, double sz);
void ndMesh3dTranslate(MESH3D * mesh3d, double tx, double ty, double tz);
void ndMesh3dRotate(MESH3D * mesh3d, double ux, double uy, double uz, double angle);
void ndMesh3dTransform(MESH3D * mesh3d, double * array_transform);


#define OBJECT3D_MESH3D								0
#define OBJECT3D_CURVE3D							1

long ndObject3dCreate(long type);
void ndObject3dSetName(long id_object3d, char * name);
void ndObject3dSetColor(long id_object3d, double r, double g, double b);
void ndObject3dSetMesh3d(long id_object3d, MESH3D * mesh3d);
void ndObject3dSetCurve3d(long id_object3d, PT3D * array_pt3d, long nb_pt3d);
void ndObject3dSetTexture(long id_object3d, unsigned char * img_rgb, long width, long height);
void ndObject3dSetWrapping(long id_object3d, double * array_tx, double * array_ty, long direction);


/*80*************************************************************************************/	
//
// Graphes : contenu
//
/*80*************************************************************************************/

void attDataGraphNew(long id_object, char * name, double x1, double y1, double x2, double y2);
void attDataGraphDelete(long id_object, char * name);
void attDataGraphSetPosition(long id_object, char * name, double x1, double y1, double x2, double y2);
void attDataGraphSetOrientation(long id_object, char * name, double angle);
void attDataGraphSetLinear(long id_object, char * name);
void attDataGraphSetPolar(long id_object, char * name);
void attDataGraphSetAutoRange(long id_object, char * name);
void attDataGraphSetXRange(long id_object, char * name, double min, double max);
void attDataGraphSetYRange(long id_object, char * name, double min, double max);
void attDataGraphSetBackgroundVisibility(long id_object, char * name, long visible);
void attDataGraphSetBackgroundColor(long id_object, char * name, double r, double g, double b);
void attDataGraphSetBorderVisibility(long id_object, char * name, long visible);
void attDataGraphSetBorderWidth(long id_object, char * name, long width);
void attDataGraphSetBorderColor(long id_object, char * name, double r, double g, double b);
void attDataGraphSetXAxisVisibility(long id_object, char * name, long visible);
void attDataGraphSetXAxisColor(long id_object, char * name, double r, double g, double b);
void attDataGraphSetXAxisWidth(long id_object, char * name, long width);
void attDataGraphSetXAxisPosition(long id_object, char * name, double position);
void attDataGraphSetYAxisVisibility(long id_object, char * name, long visible);
void attDataGraphSetYAxisColor(long id_object, char * name, double r, double g, double b);
void attDataGraphSetYAxisWidth(long id_object, char * name, long width);
void attDataGraphSetYAxisPosition(long id_object, char * name, double position);
void attDataGraphSetXGridVisibility(long id_object, char * name, long visible);
void attDataGraphSetXGridColor(long id_object, char * name, double r, double g, double b);
void attDataGraphSetYGridVisibility(long id_object, char * name, long visible);
void attDataGraphSetYGridColor(long id_object, char * name, double r, double g, double b);
void attDataGraphSetXFineGridVisibility(long id_object, char * name, long visible);
void attDataGraphSetXFineGridColor(long id_object, char * name, double r, double g, double b);
void attDataGraphSetYFineGridVisibility(long id_object, char * name, long visible);
void attDataGraphSetYFineGridColor(long id_object, char * name, double r, double g, double b);
void attDataGraphSetData(long id_object, char * name, long no, double * x, double * y, long nb);
void attDataGraphSetDataVisibility(long id_object, char * name, long no, long visible);
void attDataGraphSetDataWidth(long id_object, char * name, long no, long width);
void attDataGraphSetDataColor(long id_object, char * name, long no, double r, double g, double b);


#include <callback_prev.h>

#endif






