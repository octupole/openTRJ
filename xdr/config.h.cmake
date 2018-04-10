/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Define if building universal (internal helper macro) */
#cmakedefine AC_APPLE_UNIVERSAL_BUILD

/* Hardware and OS version for build host */
#cmakedefine BUILD_MACHINE

/* Date and time for build */
#cmakedefine BUILD_TIME

/* User doing build */
#cmakedefine BUILD_USER

/* Set to F77_FUNC(name,NAME) if Fortran used, otherwise 'name' for C. */
#cmakedefine F77_OR_C_FUNC

/* Set to F77_FUNC_(name,NAME) if Fortran used, otherwise 'name' for C. */
#cmakedefine F77_OR_C_FUNC_

/* Integer byte order is big endian. */
#cmakedefine GMX_INTEGER_BIG_ENDIAN

/* Use our own instead of system XDR libraries */
#cmakedefine GMX_INTERNAL_XDR

/* Define to 1 if the system has the type `bool'. */
#cmakedefine HAVE_BOOL

/* Define to 1 if you don't have `vprintf' but do have `_doprnt.' */
#cmakedefine HAVE_DOPRNT

/* Define to 1 if you have the `exit' function. */
#cmakedefine HAVE_EXIT

/* Define to 1 if fseeko (and presumably ftello) exists and is declared. */
#cmakedefine HAVE_FSEEKO

/* Define to 1 if you have the <inttypes.h> header file. */
#cmakedefine HAVE_INTTYPES_H

/* Define to 1 if you have the `nsl' library (-lnsl). */
#cmakedefine HAVE_LIBNSL

/* Define to 1 if you have the <memory.h> header file. */
#cmakedefine HAVE_MEMORY_H

/* Define to 1 if you have the `mkstemp' function. */
#cmakedefine HAVE_MKSTEMP

/* Define to 1 if you have the `remove' function. */
#cmakedefine HAVE_REMOVE

/* Define to 1 if you have the <rpc/xdr.h> header file. */
#cmakedefine HAVE_RPC_XDR_H

/* Define to 1 if you have the <rpc/rpc.h> header file. */
#cmakedefine HAVE_RPC_RPC_H

#cmakedefine BABA

/* Define to 1 if you have the `sqrt' function. */
#cmakedefine HAVE_SQRT

/* Define to 1 if stdbool.h conforms to C99. */
#cmakedefine HAVE_STDBOOL_H

/* Define to 1 if you have the <stdint.h> header file. */
#cmakedefine HAVE_STDINT_H

/* Define to 1 if you have the <stdlib.h> header file. */
#cmakedefine HAVE_STDLIB_H

/* Define to 1 if you have the `strcasecmp' function. */
#cmakedefine HAVE_STRCASECMP

/* Define to 1 if you have the `strdup' function. */
#cmakedefine HAVE_STRDUP

/* Define to 1 if you have the <strings.h> header file. */
#cmakedefine HAVE_STRINGS_H

/* Define to 1 if you have the <string.h> header file. */
#cmakedefine HAVE_STRING_H

/* Define to 1 if you have the <sys/stat.h> header file. */
#cmakedefine HAVE_SYS_STAT_H

/* Define to 1 if you have the <sys/types.h> header file. */
#cmakedefine HAVE_SYS_TYPES_H

/* Define to 1 if you have the <unistd.h> header file. */
#cmakedefine HAVE_UNISTD_H

/* Define to 1 if you have the `vprintf' function. */
#cmakedefine HAVE_VPRINTF

/* Define to 1 if the system has the type `_Bool'. */
#cmakedefine HAVE__BOOL

/* Define to 1 if your C compiler doesn't accept -c and -o together. */
#cmakedefine NO_MINUS_C_MINUS_O

/* Name of package */
#cmakedefine PACKAGE

/* Define to the address where bug reports for this package should be sent. */
#cmakedefine PACKAGE_BUGREPORT

/* Define to the full name of this package. */
#cmakedefine PACKAGE_NAME

/* Define to the full name and version of this package. */
#cmakedefine PACKAGE_STRING

/* Define to the one symbol short name of this package. */
#cmakedefine PACKAGE_TARNAME

/* Define to the home page for this package. */
#cmakedefine PACKAGE_URL

/* Define to the version of this package. */
#cmakedefine PACKAGE_VERSION

/* Define as the return type of signal handlers (`int' or `void'). */
#cmakedefine RETSIGTYPE

/* The size of `int', as computed by sizeof. */
#cmakedefine SIZEOF_INT ${SIZEOF_INT}

/* The size of `long int', as computed by sizeof. */
#cmakedefine SIZEOF_LONG_INT ${SIZEOF_LONG_INT}

/* The size of `long long int', as computed by sizeof. */
#cmakedefine SIZEOF_LONG_LONG_INT ${SIZEOF_LONG_LONG_INT}

/* Define to 1 if you have the ANSI C header files. */
#cmakedefine STDC_HEADERS

/* Version number of package */
#cmakedefine VERSION

/* Define if using the dmalloc debugging malloc package */
#cmakedefine WITH_DMALLOC


/* Enable large inode numbers on Mac OS X 10.5.  */
#ifndef _DARWIN_USE_64_BIT_INODE
# define _DARWIN_USE_64_BIT_INODE 1
#endif

/* Number of bits in a file offset, on hosts where this is settable. */
#cmakedefine _FILE_OFFSET_BITS

/* Define to 1 to make fseeko visible on some hosts (e.g. glibc 2.2). */
#cmakedefine _LARGEFILE_SOURCE

/* Define for large files, on AIX-style hosts. */
#cmakedefine _LARGE_FILES

/* Enable MPI compilation */
#cmakedefine __MPI

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
#cmakedefine inline
#endif

/* Define to `unsigned int' if <sys/types.h> does not define. */
#cmakedefine size_t
