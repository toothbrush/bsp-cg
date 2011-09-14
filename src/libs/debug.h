// ---- BEGIN DEBUG OUTPUT ----
#define STRINGIFY( in ) #in
#define MACROTOSTRING( in ) STRINGIFY( in )
//use the AT macro to insert the current file name and line
#define AT __FILE__ ":" MACROTOSTRING(__LINE__)
#define HERE_NOP( ... ) out ( -1 , AT, __VA_ARGS__ )
#define HERE( ... ) out( s, AT, __VA_ARGS__ )
// ---- END DEBUG OUTPUT ----

