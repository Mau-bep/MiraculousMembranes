#ifndef F77_HEADER_INCLUDED
#define F77_HEADER_INCLUDED

/* Mangling for Fortran global symbols without underscores. */
#define F77_GLOBAL(name,NAME) name##_

/* Mangling for Fortran global symbols with underscores. */
#define F77_GLOBAL_(name,NAME) name##_

#endif
