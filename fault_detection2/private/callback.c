
//#include <stdio.h>
#include <stdarg.h>

#include <callback.h>

#undef RET_LONG
#define RET_LONG		return (long)(( long(*)

#undef RET_PVOID
#define RET_PVOID		return (void *)(( void *(*)


/*80*************************************************************************************/	
//
// Entrées : déclaration
//
/*80*************************************************************************************/

long ndImageGetIn()
{
	RET_LONG(void) \
	ATT_CALL("ndImageGetIn")))();
}


long ndBlockGetIn()
{
	RET_LONG(void) \
	ATT_CALL("ndBlockGetIn")))();
}


long ndStudyGetIn()
{
	RET_LONG(void) \
	ATT_CALL("ndStudyGetIn")))();
}


long ndView3dGetIn()
{
	RET_LONG(void) \
	ATT_CALL("ndView3dGetIn")))();
}


/*80*************************************************************************************/	
//
// Sorties : déclaration, initialisation et création
//
/*80*************************************************************************************/

long ndImageNew()
{
  RET_LONG(void) \
  ATT_CALL("ndImageNew")))();
}


long ndBlockNew()
{
	RET_LONG(void) \
	ATT_CALL("ndBlockNew")))();
}


long ndGraphNew()
{
	RET_LONG(void) \
	ATT_CALL("ndGraphNew")))();
}


long ndView3dNew()
{
	RET_LONG(void) \
	ATT_CALL("ndView3dNew")))();
}


void ndImageSetTypeData(long id, long type_data)
{
	RET_VOID(long id, long type_data) \
	ATT_CALL("ndImageSetTypeData")))(id, type_data);
}


void ndImageSetSize(long id, long width, long height)
{
	RET_VOID(long id, long width, long height) \
	ATT_CALL("ndImageSetSize")))(id, width, height);
}


void ndBlockSetTypeData(long id, long type_data)
{
	RET_VOID(long id, long type_data) \
	ATT_CALL("ndBlockSetTypeData")))(id, type_data);
}


void ndBlockSetSize(long id, long width, long height, long depth)
{
	RET_VOID(long id, long width, long height, long depth) \
	ATT_CALL("ndBlockSetSize")))(id, width, height, depth);
}


void ndBlockSetFile(long id, long type_file, char * filename)
{
	RET_VOID(long id, long type_file, char * filename) \
	ATT_CALL("ndBlockSetFile")))(id, type_file, filename);
}


void ndBlockSetMemory(long id)
{
	RET_VOID(long id) \
	ATT_CALL("ndBlockSetMemory")))(id);
}


void ndGraphSetSize(long id, long width, long height)
{
	RET_VOID(long id, long width, long height) \
	ATT_CALL("ndGraphSetSize")))(id, width, height);
}


void ndImageShow(long id)
{
  RET_VOID(long id) \
	ATT_CALL("ndImageShow")))(id);
}


void ndBlockShow(long id)
{
  RET_VOID(long id) \
	ATT_CALL("ndBlockShow")))(id);
}


void ndGraphShow(long id)
{
  RET_VOID(long id) \
	ATT_CALL("ndGraphShow")))(id);
}


void ndView3dShow(long id)
{
	RET_VOID(long id) \
	ATT_CALL("ndView3dShow")))(id);
}


/*80*************************************************************************************/	
//
// Images 2D
//
/*80*************************************************************************************/

void ndImageGetSize(long id, long * width, long * height)
{
	RET_VOID(long id, long * width, long * height) \
	ATT_CALL("ndImageGetSize")))(id, width, height);
}


long ndImageGetClickNumber(long id)
{
	RET_LONG(long id) \
	ATT_CALL("ndImageGetClickNumber")))(id);
}


void ndImageGetClickArray(long id, long ** array_x, long ** array_y)
{
  RET_VOID(long id, long ** array_x, long ** array_y) \
	ATT_CALL("ndImageGetClickArray")))(id, array_x, array_y);
}


long ndImageGetCaptureNumber(long id)
{
  RET_LONG(long id) \
	ATT_CALL("ndImageGetCaptureNumber")))(id);
}


void ndImageGetCaptureArray(long id, long ** array_x, long ** array_y)
{
  RET_VOID(long id, long ** array_x, long ** array_y) \
	ATT_CALL("ndImageGetCaptureArray")))(id, array_x, array_y);
}


void * ndImageGetData(long id, long type_data)
{
  RET_PVOID(long id, long type_data) \
	ATT_CALL("ndImageGetData")))(id, type_data);
}


void * ndImageGetPtrData(long id)
{
  RET_PVOID(long id) \
	ATT_CALL("ndImageGetPtrData")))(id);
}


long ndImageGetTypeData(long id)
{
  RET_LONG(long id) \
	ATT_CALL("ndImageGetTypeData")))(id);
}


void ndImageSetClass(long id, long type_class)
{
  RET_VOID(long id, long type_class) \
	ATT_CALL("ndImageSetClass")))(id, type_class);
}


void ndImageAddObject2d(long id, long id_object2d)
{
  RET_VOID(long id, long id_object2d) \
	ATT_CALL("ndImageAddObject2d")))(id, id_object2d);
}


void ndImageDeleteObject2d(long id, long id_object2d)
{
  RET_VOID(long id, long id_object2d) \
	ATT_CALL("ndImageDeleteObject2d")))(id, id_object2d);
}


void ndImageDeleteObject2dByName(long id, char * name)
{
  RET_VOID(long id, char * name) \
	ATT_CALL("ndImageDeleteObject2dByName")))(id, name);
}


long ndImageAddLine(long id, double x1, double y1, double x2, double y2)
{
  RET_LONG(long id, double x1, double y1, double x2, double y2) \
	ATT_CALL("ndImageAddLine")))(id, x1, y1, x2, y2);
}

long ndImageAddSquare(long id, double x, double y, double c)
{
  RET_LONG(long id, double x, double y, double c) \
	ATT_CALL("ndImageAddSquare")))(id, x, y, x);
}


long ndImageSaveBitmap(long id, char * filename)
{
  RET_LONG(long id, char * filename) \
	ATT_CALL("ndImageSaveBitmap")))(id, filename);
}


/*80*************************************************************************************/	
//
// Images 3D
//
/*80*************************************************************************************/

void ndBlockGetSize(long id, long * width, long * height, long * depth)
{
  RET_VOID(long id, long * width, long * height, long * depth) \
	ATT_CALL("ndBlockGetSize")))(id, width, height, depth);
}


void * ndBlockGetFront(long id, long no, long type_data)
{
  RET_PVOID(long id, long no, long type_data) \
	ATT_CALL("ndBlockGetFront")))(id, no, type_data);
}


void * ndBlockGetPartialFront(long id, long no, long type_data, long x, long y, long dx, long dy)
{
  RET_PVOID(long id, long no, long type_data, long x, long y, long dx, long dy) \
	ATT_CALL("ndBlockGetPartialFront")))(id, no, type_data, x, y, dx, dy);
}


void * ndBlockGetRight(long id, long no, long type_data)
{
  RET_PVOID(long id, long no, long type_data) \
	ATT_CALL("ndBlockGetRight")))(id, no, type_data);
}


void * ndBlockGetPartialRight(long id, long no, long type_data, long x, long y, long dx, long dy)
{
  RET_PVOID(long id, long no, long type_data, long x, long y, long dx, long dy) \
	ATT_CALL("ndBlockGetPartialRight")))(id, no, type_data, x, y, dx, dy);
}


void * ndBlockGetTop(long id, long no, long type_data)
{
  RET_PVOID(long id, long no, long type_data) \
	ATT_CALL("ndBlockGetTop")))(id, no, type_data);
}


void * ndBlockGetPartialTop(long id, long no, long type_data, long x, long y, long dx, long dy)
{
	RET_PVOID(long id, long no, long type_data, long x, long y, long dx, long dy) \
	ATT_CALL("ndBlockGetPartialTop")))(id, no, type_data, x, y, dx, dy);
}


void * ndBlockGetImageFromTopPatch(long id, long no, long type_data, long * patch_dx, long * patch_dy, long nb, long * width, long * height)
{
  RET_PVOID(long id, long no, long type_data, long * patch_dx, long * patch_dy, long nb, long * width, long * height) \
	ATT_CALL("ndBlockGetImageFromTopPatch")))(id, no, type_data, patch_dx, patch_dy, nb, width, height);
}


void * ndBlockGetImageSONE(long id, long no, long type_data, long * width, long * height)
{
  RET_PVOID(long id, long no, long type_data, long * width, long * height) \
	ATT_CALL("ndBlockGetImageSONE")))(id, no, type_data, width, height);
}


void * ndBlockGetImageNOSE(long id, long no, long type_data, long * width, long * height)
{
  RET_PVOID(long id, long no, long type_data, long * width, long * height) \
	ATT_CALL("ndBlockGetImageNOSE")))(id, no, type_data, width, height);
}


void ndBlockSetUpdateFile(long id, long type_update)
{
	RET_VOID(long id, long type_update) \
	ATT_CALL("ndBlockSetUpdateFile")))(id, type_update);
}


void ndBlockSetPartialFront(long id, long no, long x, long y, void * data, long w, long h)
{
	RET_VOID(long id, long no, long x, long y, void * data, long w, long h) \
	ATT_CALL("ndBlockSetPartialFront")))(id, no, x, y, data, w, h);
}


void ndBlockSetFront(long id, long no, void * data)
{
	RET_VOID(long id, long no, void * data) \
	ATT_CALL("ndBlockSetFront")))(id, no, data);
}


void ndBlockSetFrontByType(long id, long no, void * data, long type_data)
{
  RET_VOID(long id, long no, void * data, long type_data) \
	ATT_CALL("ndBlockSetFrontByType")))(id, no, data, type_data);
}


void ndBlockSetPartialFrontByType(long id, long no, long x, long y, void * data, long w, long h, long type_data)
{
  RET_VOID(long id, long no, long x, long y, void * data, long w, long h, long type_data) \
	ATT_CALL("ndBlockSetPartialFrontByType")))(id, no, x, y, data, w, h, type_data);
}


void ndBlockSetPartialRight(long id, long no, long x, long y, void * data, long w, long h)
{
  RET_VOID(long id, long no, long x, long y, void * data, long w, long h) \
	ATT_CALL("ndBlockSetPartialRight")))(id, no, x, y, data, w, h);
}


void ndBlockSetRight(long id, long no, void * data)
{
  RET_VOID(long id, long no, void * data) \
	ATT_CALL("ndBlockSetRight")))(id, no, data);
}


void ndBlockSetRightByType(long id, long no, void * data, long type_data)
{
  RET_VOID(long id, long no, void * data, long type_data) \
	ATT_CALL("ndBlockSetRightByType")))(id, no, data, type_data);
}


void ndBlockSetPartialRightByType(long id, long no, long x, long y, void * data, long w, long h, long type_data)
{
  RET_VOID(long id, long no, long x, long y, void * data, long w, long h, long type_data) \
	ATT_CALL("ndBlockSetPartialRightByType")))(id, no, x, y, data, w, h, type_data);
}


void ndBlockSetPartialTop(long id, long no, long x, long y, void * data, long w, long h)
{
  RET_VOID(long id, long no, long x, long y, void * data, long w, long h) \
	ATT_CALL("ndBlockSetPartialTop")))(id, no, x, y, data, w, h);
}

void ndBlockSetTop(long id, long no, void * data)
{
	RET_VOID(long id, long no, void * data) \
	ATT_CALL("ndBlockSetTop")))(id, no, data);
}


void ndBlockSetTopByType(long id, long no, void * data, long type_data)
{
  RET_VOID(long id, long no, void * data, long type_data) \
	ATT_CALL("ndBlockSetTopByType")))(id, no, data, type_data);
}


void ndBlockSetPartialTopByType(long id, long no, long x, long y, void * data, long w, long h, long type_data)
{
  RET_VOID(long id, long no, long x, long y, void * data, long w, long h, long type_data) \
	ATT_CALL("ndBlockSetPartialTopByType")))(id, no, x, y, data, w, h, type_data);
}


void ndBlockSetImageSONE(long id, long no, void * data)
{
  RET_VOID(long id, long no, void * data) \
	ATT_CALL("ndBlockSetImageSONE")))(id, no, data);
}


void ndBlockSetImageSONEByType(long id, long no, void * data, long type_data)
{
  RET_VOID(long id, long no, void * data, long type_data) \
	ATT_CALL("ndBlockSetImageSONEByType")))(id, no, data, type_data);
}


void ndBlockSetImageNOSE(long id, long no, void * data)
{
  RET_VOID(long id, long no, void * data) \
	ATT_CALL("ndBlockSetImageNOSE")))(id, no, data);
}


void ndBlockSetImageNOSEByType(long id, long no, void * data, long type_data)
{
  RET_VOID(long id, long no, void * data, long type_data) \
	ATT_CALL("ndBlockSetImageNOSEByType")))(id, no, data, type_data);
}


long ndBlockGetClickNumber(long id)
{
  RET_LONG(long id) \
	ATT_CALL("ndBlockGetClickNumber")))(id);
}


void ndBlockGetClickArray(long id, long ** array_x, long ** array_y, long ** array_z)
{
  RET_VOID(long id, long ** array_x, long ** array_y, long ** array_z) \
	ATT_CALL("ndBlockGetClickArray")))(id, array_x, array_y, array_z);
}


long ndBlockGetCaptureNumber(long id)
{
  RET_LONG(long id) \
	ATT_CALL("ndBlockGetCaptureNumber")))(id);
}


void ndBlockGetCaptureArray(long id, long ** array_x, long ** array_y, long ** array_z)
{
  RET_VOID(long id, long ** array_x, long ** array_y, long ** array_z) \
	ATT_CALL("ndBlockGetCaptureArray")))(id, array_x, array_y, array_z);
}


long ndBlockGetCursor(long id, long * x, long * y, long * z)
{
  RET_LONG(long id, long * x, long * y, long * z) \
	ATT_CALL("ndBlockGetCursor")))(id, x, y, z);
}


void * ndBlockGetPtrData(long id)
{
  RET_PVOID(long id) \
	ATT_CALL("ndBlockGetPtrData")))(id);
}


long ndBlockGetTypeData(long id)
{
  RET_LONG(long id) \
	ATT_CALL("ndBlockGetTypeData")))(id);
}


void ndBlockSetClass(long id, long type_class)
{
  RET_VOID(long id, long type_class) \
	ATT_CALL("ndBlockSetClass")))(id, type_class);
}


void ndBlockCopyToFile(long id, char * filename)
{
	RET_VOID(long id, char * filename) \
	ATT_CALL("ndBlockCopyToFile")))(id, filename);
}


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

void ndView3dGetSize(long id, long * width, long * height, long * depth)
{
	RET_VOID(long id, long * width, long * height, long * depth) \
	ATT_CALL("ndView3dGetSize")))(id, width, height, depth);
}


void ndView3dSetSize(long id, long width, long height, long depth)
{
  RET_VOID(long id, long width, long height, long depth) \
	ATT_CALL("ndView3dSetSize")))(id, width, height, depth);
}


void ndView3dAddObject3d(long id, long id_object3d)
{
  RET_VOID(long id, long id_object3d) \
	ATT_CALL("ndView3dAddObject3d")))(id, id_object3d);
}


void ndView3dDeleteObject3d(long id, long id_object3d)
{
  RET_VOID(long id, long id_object3d) \
	ATT_CALL("ndView3dDeleteObject3d")))(id, id_object3d);
}


void ndView3dDeleteObject3dByName(long id, char * name)
{
  RET_VOID(long id, char * name) \
	ATT_CALL("ndView3dDeleteObject3dByName")))(id, name);
}

void ndObject3dSetTransparency(long id, double alpha)
{
  RET_VOID(long id_object3d, double alpha) \
	ATT_CALL("ndObject3dSetTransparency")))(id, alpha);
}


/*80*************************************************************************************/	
//
// Vidéos
//
/*80*************************************************************************************/

long ndVideoGetId()
{
	RET_LONG(void) \
	ATT_CALL("ndVideoGetId")))();
}


void ndVideoSetFormatIn(long id, long format)
{
  RET_VOID(long id, long format) \
	ATT_CALL("ndVideoSetFormatIn")))(id, format);
}


void ndVideoSetFormatOut(long id, long format)
{
  RET_VOID(long id, long format) \
	ATT_CALL("ndVideoSetFormatOut")))(id, format);
}


void * ndVideoGetDataIn(long id)
{
	RET_PVOID(long id) \
	ATT_CALL("ndVideoGetDataIn")))(id);
}


void * ndVideoGetDataOut(long id)
{
  RET_PVOID(long id) \
	ATT_CALL("ndVideoGetDataOut")))(id);
}


long ndVideoGetWidth(long id)
{
  RET_LONG(long id) \
	ATT_CALL("ndVideoGetWidth")))(id);
}


long ndVideoGetHeight(long id)
{
  RET_LONG(long id) \
	ATT_CALL("ndVideoGetHeight")))(id);
}


void * ndVideoGetUserData(long id)
{
	RET_PVOID(long id) \
	ATT_CALL("ndVideoGetUserData")))(id);
}


void ndVideoSetUserData(long id, void * user_data)
{
  RET_VOID(long id, void * user_data) \
	ATT_CALL("ndVideoSetUserData")))(id, user_data);
}


void ndVideoStart(long id, ndfunc videoOn, ndfunc videoOff)
{
  RET_VOID(long id, ndfunc VideoOn, ndfunc VideoOff) \
	ATT_CALL("ndVideoStart")))(id, videoOn, videoOff);
}


long ndVideoGetDuration()
{
  RET_LONG(void) \
	ATT_CALL("ndVideoGetDuration")))();
}


long ndVideoGetCurrentDuration()
{
  RET_LONG(void) \
	ATT_CALL("ndVideoGetCurrentDuration")))();
}


long ndVideoGetCurrentTime()
{
  RET_LONG(void) \
	ATT_CALL("ndVideoGetCurrentTime")))();
}


long ndVideoGetNumberOfFrames(long id)
{
  RET_LONG(long id) \
	ATT_CALL("ndVideoGetNumberOfFrames")))(id);
}


long ndVideoGetNoCurrentFrame(long id)
{
  RET_LONG(long id) \
	ATT_CALL("ndVideoGetNoCurrentFrame")))(id);
}


void ndVideoMoveTo(long id, long no)
{
  RET_VOID(long id, long no) \
	ATT_CALL("ndVideoMoveTo")))(id, no);
}


void ndVideoRefreshOut(long id)
{
  RET_VOID(long id) \
	ATT_CALL("ndVideoRefreshOut")))(id);
}



void * ndVideoFileCreate(char * filename)
{
	RET_PVOID(char * filename) \
	ATT_CALL("ndVideoFileCreate")))(filename);
}


void ndVideoFileSet(void * avi, long width, long height, double frequency)
{
	RET_VOID(void * avi, long width, long height, double frequency) \
	ATT_CALL("ndVideoFileSet")))(avi, width, height, frequency);
}


void ndVideoFileAddFrame(void * avi, unsigned char * img_rgba, long width, long height)
{
	RET_VOID(void * avi, unsigned char * img_rgba, long width, long height) \
	ATT_CALL("ndVideoFileAddFrame")))(avi, img_rgba, width, height);
}


void * ndVideoFileClose(void * avi)
{
	RET_PVOID(void * avi) \
	ATT_CALL("ndVideoFileClose")))(avi);
}


/*80*************************************************************************************/	
//
// Fenêtres : images 2D, images 3D, vues 3D, vidéos, 
//
/*80*************************************************************************************/

void ndWindowRedraw(long id_object)
{
  RET_VOID(long id_object) \
	ATT_CALL("ndWindowRedraw")))(id_object);
}


void ndWindowClose(long id_object)
{
	RET_VOID(long hwnd) \
	ATT_CALL("attWindowClose")))(id_object);
}


void ndWindowSetName(long id_object, const char * name)
{
  RET_VOID(long id_object, const char * name) \
	ATT_CALL("ndWindowSetName")))(id_object, name);
}


void ndWindowUserDataNew(long id_object, char * name)
{
  RET_VOID(long id_object, char * name) \
	ATT_CALL("ndWindowUserDataNew")))(id_object, name);
}


void ndWindowUserDataDelete(long id_object, char * name)
{
  RET_VOID(long id_object, char * name) \
	ATT_CALL("ndWindowUserDataDelete")))(id_object, name);
}


void ndWindowUserDataSet(long id_object, char * name, void * ptr)
{
	RET_VOID(long id_object, char * name, void * ptr) \
	ATT_CALL("ndWindowUserDataSet")))(id_object, name, ptr);
}


void * ndWindowUserDataGet(long id_object, char * name)
{
  RET_PVOID(long id_object, char * name) \
	ATT_CALL("ndWindowUserDataGet")))(id_object, name);
}


void ndWindowCallbackAdd(long id_object, char * name, ndfunc_function_callback function_callback, void * data_callback)
{
  RET_VOID(long id_object, char * name, ndfunc_function_callback function_callback, void * data_callback) \
	ATT_CALL("ndWindowCallbackAdd")))(id_object, name, function_callback, data_callback);
}


void ndWindowCallbackSub(long id_object, char * name, ndfunc_function_callback function_callback)
{
  RET_VOID(long id_object, char * name, ndfunc_function_callback function_callback) \
	ATT_CALL("ndWindowCallbackSub")))(id_object, name, function_callback);
}


long ndWindowCallbackEventGet(long id_object, void * data_event, char * name_event, void * data_to_get)
{
  RET_LONG(long id_object, void * data_event, char * name_event, void * data_to_get) \
	ATT_CALL("ndWindowCallbackEventGet")))(id_object, data_event, name_event, data_to_get);
}


long ndWindowCallbackEventSet(long id_object, void * data_event, char * name_event, void * data_to_set)
{
  RET_LONG(long id_object, void * data_event, char * name_event, void * data_to_set) \
	ATT_CALL("ndWindowCallbackEventSet")))(id_object, data_event, name_event, data_to_set);
}


void * ndWindowGetHwnd(long id_object)
{
	RET_PVOID(long id_object) \
	ATT_CALL("attDataGetHwnd")))(id_object);
}


/*80*************************************************************************************/	
//
// Divers
//
/*80*************************************************************************************/

void ndProgressSetMinMax(long mini, long maxi)
{
  RET_VOID(long mini, long maxi) \
	ATT_CALL("ndProgressSetMinMax")))(mini, maxi);
}


void ndProgressSetPosition(long position)
{
  RET_VOID(long position) \
	ATT_CALL("ndProgressSetPosition")))(position);
}


void ndProgressReset()
{
  RET_VOID() \
	ATT_CALL("ndProgressReset")))();
}


long ndPopupImage(void * data, long w, long h, long type_data)
{
  RET_LONG(void * data, long w, long h, long type_data) \
	ATT_CALL("ndPopupImage")))(data, w, h, type_data);
}


long ndPopupHistogram(void * data, long length, long type_data)
{
  RET_LONG(void * data, long length, long type_data) \
	ATT_CALL("ndPopupImage")))(data, length, type_data);
}


void ndFree(void * ptr)
{
  RET_VOID(void * ptr) \
	ATT_CALL("ndFree")))(ptr);
}


void ndMsg(char * message)
{
  RET_VOID(char * message) \
	ATT_CALL("ndMsg")))(message);
}


long ndMenuCreate(char * title)
{
  RET_LONG(char * title) \
	ATT_CALL("ndMenuCreate")))(title);
}


void ndMenuAddItem(long menu, char * title, ndfunc function)
{
	RET_VOID(long menu, char * title, ndfunc function) \
	ATT_CALL("ndMenuAddItem")))(menu, title, function);
}


void ndMenuAddSeparator(long menu)
{
  RET_VOID(long menu) \
	ATT_CALL("ndMenuAddSeparator")))(menu);
}


void ndMenuAddSubMenu(long menu, long submenu)
{
  RET_VOID(long menu, long submenu) \
	ATT_CALL("ndMenuAddSubMenu")))(menu, submenu);
}


long ndBoxParameterDialog(void * data)
{
  RET_LONG(void * data) \
	ATT_CALL("ndBoxParameterDialog")))(data);
}


void ndMsgBox(char * title, char * message)
{
  RET_VOID(char * title, char * message) \
	ATT_CALL("attMsgBox")))(title, message);
}


long ndMsgBoxYNC(char * title, char * message)
{
  RET_LONG(char * title, char * message) \
	ATT_CALL("attMsgBoxYNC")))(title, message);
}


void ndLibraryAddMenu(long id_library, char * insert, long menu)
{
  RET_VOID(long id_library, char * insert, long menu) \
	ATT_CALL("ndLibraryAddMenu")))(id_library, insert, menu);
}


void ndLibraryAddMenuPopup(long id_library, long id_object, long x, long y, long menu)
{
	RET_VOID(long id_library, long id_object, long x, long y, long menu) \
	ATT_CALL("ndLibraryAddMenuPopup")))(id_library, id_object, x, y, menu);
}

//long ndLibraryFalseLoad(long id_library, char * name, attfunc functionFake);
long ndLibraryFalseLoad(long id_library, char * name, ndfunc functionFalse)
{
	RET_LONG(long id_library, char * name, ndfunc functionFalse) \
	ATT_CALL("ndLibraryFalseLoad")))(id_library, name, functionFalse);
	return 0;
}


// Répertoire
//char * ndDirectoryScanGetFiles(char * dirname) ;
//char * ndDirectoryScanGetDirectories(char * dirname);
char * ndDirectoryScanGetFiles(char * dirname)
{
	RET_PCHAR(char * dirname) \
	ATT_CALL("ndDirectoryScanGetFiles")))(dirname);
	return NULL;
}


char * ndDirectoryScanGetDirectories(char * dirname)
{
	RET_PCHAR(char * dirname) \
	ATT_CALL("ndDirectoryScanGetDirectories")))(dirname);
	return NULL;
}


void ndLoadFromFile(char * filename)
{
	RET_VOID(char * filename) \
	ATT_CALL("ndLoadFromFile")))(filename);
}

















/*80*************************************************************************************/	
//
// Objets 2D
//
/*80*************************************************************************************/

long ndObject2dCreate(long type)
{
  RET_LONG(long type) \
	ATT_CALL("ndObject2dCreate")))(type);
}


void ndObject2dSetName(long id_object2d, char * name)
{
  RET_VOID(long id_object2d, char * name) \
	ATT_CALL("ndObject2dSetName")))(id_object2d, name);
}


void ndObject2dSetColor(long id_object2d, double r, double g, double b)
{
	RET_VOID(long id_object2d, double r, double g, double b) \
	ATT_CALL("ndObject2dSetColor")))(id_object2d, r, g, b);
}


void ndObject2dSetThickness(long id_object2d, long thickness)
{
  RET_VOID(long id_object2d, long thickness) \
	ATT_CALL("ndObject2dSetThickness")))(id_object2d, thickness);
}


void ndObject2dSetArrayPoints(long id_object2d, PT2D * array_pt2d, long nb_pt2d)
{
  RET_VOID(long id_object2d, PT2D * array_pt2d, long nb_pt2d) \
	ATT_CALL("ndObject2dSetArrayPoints")))(id_object2d, array_pt2d, nb_pt2d);
}


void ndObject2dTranslate(long id_object2d, double sx, double sy)
{
  RET_VOID(long id_object2d, double sx, double sy) \
	ATT_CALL("ndObject2dTranslate")))(id_object2d, sx, sy);
}


void ndObject2dScale(long id_object2d, double sx, double sy)
{
  RET_VOID(long id_object2d, double sx, double sy) \
	ATT_CALL("ndObject2dScale")))(id_object2d, sx, sy);
}


/*80*************************************************************************************/	
//
// Objets 3D
//
/*80*************************************************************************************/

void * ndMesh3dCube()
{
  RET_PVOID(void) \
	ATT_CALL("ndMesh3dCube")))();
}


void * ndMesh3dSphere(long precision)
{
	RET_PVOID(long precision) \
	ATT_CALL("ndMesh3dSphere")))(precision);
}


void * ndMesh3dCylinder(long precision)
{
  RET_PVOID(long precision) \
	ATT_CALL("ndMesh3dCylinder")))(precision);
}


void * ndMesh3dCone(long precision)
{
  RET_PVOID(long precision) \
	ATT_CALL("ndMesh3dCone")))(precision);
}


void * ndMesh3dArrow(PT3D * pt_begin, PT3D * pt_end, double radius)
{
  RET_PVOID(PT3D * pt_begin, PT3D * pt_end, double radius) \
	ATT_CALL("ndMesh3dArrow")))(pt_begin, pt_end, radius);
}


void * ndMesh3dLine(PT3D * array_pt, long nb_pt, double radius, long precision)
{
  RET_PVOID(PT3D * array_pt, long nb_pt, double radius, long precision) \
	ATT_CALL("ndMesh3dLine")))(array_pt, nb_pt, radius, precision);
}


void * ndMesh3dReadFromFile(char * filename)
{
	RET_PVOID(char * filename) \
	ATT_CALL("ndMesh3dReadFromFile")))(filename);
}


void * ndMesh3dAlloc(long nb_pt, long nb_triangle)
{
  RET_PVOID(long nb_pt, long nb_triangle) \
	ATT_CALL("ndMesh3dAlloc")))(nb_pt, nb_triangle);
}


void * ndMesh3dFree(MESH3D * mesh3d)
{
  RET_PVOID(MESH3D * mesh3d) \
	ATT_CALL("ndMesh3dFree")))(mesh3d);
}


void ndMesh3dScale(MESH3D * mesh3d, double sx, double sy, double sz)
{
  RET_VOID(MESH3D * mesh3d, double sx, double sy, double sz) \
	ATT_CALL("ndMesh3dScale")))(mesh3d, sx, sy, sz);
}


void ndMesh3dTranslate(MESH3D * mesh3d, double tx, double ty, double tz)
{
  RET_VOID(MESH3D * mesh3d, double tx, double ty, double tz) \
	ATT_CALL("ndMesh3dTranslate")))(mesh3d, tx, ty, tz);
}


void ndMesh3dRotate(MESH3D * mesh3d, double ux, double uy, double uz, double angle)
{
  RET_VOID(MESH3D * mesh3d, double ux, double uy, double uz, double angle) \
	ATT_CALL("ndMesh3dRotate")))(mesh3d, ux, uy, uz, angle);
}


void ndMesh3dTransform(MESH3D * mesh3d, double * array_transform)
{
  RET_VOID(MESH3D * mesh3d, double * array_transform) \
	ATT_CALL("ndMesh3dTransform")))(mesh3d, array_transform);
}


long ndObject3dCreate(long type)
{
  RET_LONG(long type) \
	ATT_CALL("ndObject3dCreate")))(type);
}


void ndObject3dSetName(long id_object3d, char * name)
{
  RET_VOID(long id_object3d, char * name) \
	ATT_CALL("ndObject3dSetName")))(id_object3d, name);
}


void ndObject3dSetColor(long id_object3d, double r, double g, double b)
{
  RET_VOID(long id_object3d, double r, double g, double b) \
	ATT_CALL("ndObject3dSetColor")))(id_object3d, r, g, b);
}


void ndObject3dSetMesh3d(long id_object3d, MESH3D * mesh3d)
{
  RET_VOID(long id_object3d, MESH3D * mesh3d) \
	ATT_CALL("ndObject3dSetMesh3d")))(id_object3d, mesh3d);
}


void ndObject3dSetCurve3d(long id_object3d, PT3D * array_pt3d, long nb_pt3d)
{
  RET_VOID(long id_object3d, PT3D * array_pt3d, long nb_pt3d) \
	ATT_CALL("ndObject3dSetCurve3d")))(id_object3d, array_pt3d, nb_pt3d);
}


void ndObject3dSetTexture(long id_object3d, unsigned char * img_rgb, long width, long height)
{
  RET_VOID(long id_object3d, unsigned char * img_rgb, long width, long height) \
	ATT_CALL("ndObject3dSetTexture")))(id_object3d, img_rgb, width, height);
}


void ndObject3dSetWrapping(long id_object3d, double * array_tx, double * array_ty, long direction)
{
  RET_VOID(long id_object3d, double * array_tx, double * array_ty, long direction) \
	ATT_CALL("ndObject3dSetWrapping")))(id_object3d, array_tx, array_ty, direction);
}












/*80*************************************************************************************/	
//
// Graphes : contenu
//
/*80*************************************************************************************/

void attDataGraphNew(long id_object, char * name, double x1, double y1, double x2, double y2)
{
  RET_VOID(long id_object, char * name, double x1, double y1, double x2, double y2) \
	ATT_CALL("attDataGraphNew")))(id_object, name, x1, y1, x2, y2);
}


void attDataGraphDelete(long id_object, char * name)
{
  RET_VOID(long id_object, char * name) \
	ATT_CALL("attDataGraphDelete")))(id_object, name);
}


void attDataGraphSetPosition(long id_object, char * name, double x1, double y1, double x2, double y2)
{
  RET_VOID(long id_object, char * name, double x1, double y1, double x2, double y2) \
	ATT_CALL("attDataGraphSetPosition")))(id_object, name, x1, y1, x2, y2);
}


void attDataGraphSetOrientation(long id_object, char * name, double angle)
{
  RET_VOID(long id_object, char * name, double angle) \
	ATT_CALL("attDataGraphSetOrientation")))(id_object, name, angle);
}


void attDataGraphSetLinear(long id_object, char * name)
{
  RET_VOID(long id_object, char * name) \
	ATT_CALL("attDataGraphSetLinear")))(id_object, name);
}


void attDataGraphSetPolar(long id_object, char * name)
{
  RET_VOID(long id_object, char * name) \
	ATT_CALL("attDataGraphSetPolar")))(id_object, name);
}


void attDataGraphSetAutoRange(long id_object, char * name)
{
  RET_VOID(long id_object, char * name) \
	ATT_CALL("attDataGraphSetAutoRange")))(id_object, name);
}


void attDataGraphSetXRange(long id_object, char * name, double min, double max)
{
  RET_VOID(long id_object, char * name, double min, double max) \
	ATT_CALL("attDataGraphSetXRange")))(id_object, name, min, max);
}


void attDataGraphSetYRange(long id_object, char * name, double min, double max)
{
  RET_VOID(long id_object, char * name, double min, double max) \
	ATT_CALL("attDataGraphSetYRange")))(id_object, name, min, max);
}


void attDataGraphSetBackgroundVisibility(long id_object, char * name, long visible)
{
  RET_VOID(long id_object, char * name, long visible) \
	ATT_CALL("attDataGraphSetBackgroundVisibility")))(id_object, name, visible);
}


void attDataGraphSetBackgroundColor(long id_object, char * name, double r, double g, double b)
{
  RET_VOID(long id_object, char * name, double r, double g, double b) \
	ATT_CALL("attDataGraphSetBackgroundColor")))(id_object, name, r, g, b);
}


void attDataGraphSetBorderVisibility(long id_object, char * name, long visible)
{
  RET_VOID(long id_object, char * name, long visible) \
	ATT_CALL("attDataGraphSetBorderVisibility")))(id_object, name, visible);
}


void attDataGraphSetBorderWidth(long id_object, char * name, long width)
{
  RET_VOID(long id_object, char * name, long width) \
	ATT_CALL("attDataGraphSetBorderWidth")))(id_object, name, width);
}


void attDataGraphSetBorderColor(long id_object, char * name, double r, double g, double b)
{
  RET_VOID(long id_object, char * name, double r, double g, double b) \
	ATT_CALL("attDataGraphSetBorderColor")))(id_object, name, r, g, b);
}


void attDataGraphSetXAxisVisibility(long id_object, char * name, long visible)
{
  RET_VOID(long id_object, char * name, long visible) \
	ATT_CALL("attDataGraphSetXAxisVisibility")))(id_object, name, visible);
}


void attDataGraphSetXAxisColor(long id_object, char * name, double r, double g, double b)
{
  RET_VOID(long id_object, char * name, double r, double g, double b) \
	ATT_CALL("attDataGraphSetXAxisColor")))(id_object, name, r, g, b);
}


void attDataGraphSetXAxisWidth(long id_object, char * name, long width)
{
  RET_VOID(long id_object, char * name, long width) \
	ATT_CALL("attDataGraphSetXAxisWidth")))(id_object, name, width);
}


void attDataGraphSetXAxisPosition(long id_object, char * name, double position)
{
  RET_VOID(long id_object, char * name, double position) \
  ATT_CALL("attDataGraphSetXAxisPosition")))(id_object, name, position);
}


void attDataGraphSetYAxisVisibility(long id_object, char * name, long visible)
{
  RET_VOID(long id_object, char * name, long visible) \
	ATT_CALL("attDataGraphSetYAxisVisibility")))(id_object, name, visible);
}


void attDataGraphSetYAxisColor(long id_object, char * name, double r, double g, double b)
{
  RET_VOID(long id_object, char * name, double r, double g, double b) \
	ATT_CALL("attDataGraphSetYAxisColor")))(id_object, name, r, g, b);
}


void attDataGraphSetYAxisWidth(long id_object, char * name, long width)
{
  RET_VOID(long id_object, char * name, long width) \
	ATT_CALL("attDataGraphSetYAxisWidth")))(id_object, name, width);
}


void attDataGraphSetYAxisPosition(long id_object, char * name, double position)
{
  RET_VOID(long id_object, char * name, double position) \
	ATT_CALL("attDataGraphSetYAxisPosition")))(id_object, name, position);
}


void attDataGraphSetXGridVisibility(long id_object, char * name, long visible)
{
  RET_VOID(long id_object, char * name, long visible) \
	ATT_CALL("attDataGraphSetXGridVisibility")))(id_object, name, visible);
}


void attDataGraphSetXGridColor(long id_object, char * name, double r, double g, double b)
{
  RET_VOID(long id_object, char * name, double r, double g, double b) \
	ATT_CALL("attDataGraphSetXGridColor")))(id_object, name, r, g, b);
}


void attDataGraphSetYGridVisibility(long id_object, char * name, long visible)
{
  RET_VOID(long id_object, char * name, long visible) \
	ATT_CALL("attDataGraphSetYGridVisibility")))(id_object, name, visible);
}


void attDataGraphSetYGridColor(long id_object, char * name, double r, double g, double b)
{
  RET_VOID(long id_object, char * name, double r, double g, double b) \
	ATT_CALL("attDataGraphSetYGridColor")))(id_object, name, r, g, b);
}


void attDataGraphSetXFineGridVisibility(long id_object, char * name, long visible)
{
  RET_VOID(long id_object, char * name, long visible) \
	ATT_CALL("attDataGraphSetXFineGridVisibility")))(id_object, name, visible);
}


void attDataGraphSetXFineGridColor(long id_object, char * name, double r, double g, double b)
{
	RET_VOID(long id_object, char * name, double r, double g, double b) \
	ATT_CALL("attDataGraphSetXFineGridColor")))(id_object, name, r, g, b);
}


void attDataGraphSetYFineGridVisibility(long id_object, char * name, long visible)
{
  RET_VOID(long id_object, char * name, long visible) \
	ATT_CALL("attDataGraphSetYFineGridVisibility")))(id_object, name, visible);
}


void attDataGraphSetYFineGridColor(long id_object, char * name, double r, double g, double b)
{
  RET_VOID(long id_object, char * name, double r, double g, double b) \
	ATT_CALL("attDataGraphSetYFineGridColor")))(id_object, name, r, g, b);
}


void attDataGraphSetData(long id_object, char * name, long no, double * x, double * y, long nb)
{
  RET_VOID(long id_object, char * name, long no, double * x, double * y, long nb) \
	ATT_CALL("attDataGraphSetData")))(id_object, name, no, x, y, nb);
}


void attDataGraphSetDataVisibility(long id_object, char * name, long no, long visible)
{
  RET_VOID(long id_object, char * name, long no, long visible) \
	ATT_CALL("attDataGraphSetDataVisibility")))(id_object, name, no, visible);
}


void attDataGraphSetDataWidth(long id_object, char * name, long no, long width)
{
  RET_VOID(long id_object, char * name, long no, long width) \
	ATT_CALL("attDataGraphSetDataWidth")))(id_object, name, no, width);
}


void attDataGraphSetDataColor(long id_object, char * name, long no, double r, double g, double b)
{
	RET_VOID(long id_object, char * name, long no, double r, double g, double b) \
	ATT_CALL("attDataGraphSetDataColor")))(id_object, name, no, r, g, b);
}







