
/*80*************************************************************************************/
/*80*************************************************************************************/
/*80********************************** NOT TO MODIFY ************************************/
/*80*************************************************************************************/
/*80*************************************************************************************/

#include <windows.h>

#define LIBRARY_MAIN
#include <callback.h>


int WINAPI DllMain(HINSTANCE hInstance, DWORD fdwReason, PVOID pvReserved)
{	
	return TRUE;
}


EXPORT_LIB void GetLibraryVersion(char * version)
{
	strcpy(version, "1.0.11");
}


EXPORT_LIB SetCallbackFunction(attfunc_functionCallback callback)
{
	functionCallback = callback;
}
