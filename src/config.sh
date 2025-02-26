#! /bin/sh -f

test_include_folders="/ /usr/include /usr/local/include /sw/include"
test_lib_folders="/ /usr/lib /usr/local/lib /sw/lib"

mkdir -p "../bin"

cat > tst.c  <<EOF
#include <fitsio.h>
int main(int argc,char** argv) {
return 0;
}
EOF

FITSIO_INCLUDE_FLAG=""
for i in $test_include_folders; do
    if test $i == "/"; then
	inc_flag=""
    else
	inc_flag=-I${i}
    fi

    status=`gcc ${inc_flag} tst.c -o tst.exe ; echo $?`

    if test $status -eq 0; then
	    FITSIO_INC_FLAG=$inc_flag
	    break
    fi
done
rm -f tst.c tst.exe


cat >tst.c  <<EOF
#include <fitsio.h>
int main(int argc,char** argv) {
float version;
fits_get_version(&version);
fprintf(stdout,"%f\n",version);
return 0;
}
EOF

FITSIO_LIB_FLAG="?"
for i in $test_lib_folders;do 
    if test $i == "/"; then 
	lib_flag=""
    else
	lib_flag=-L${i}
    fi

    status=`gcc ${FITSIO_INCLUDE_FLAG} ${lib_flag} tst.c -lcfitsio -o tst.exe ; echo $?`

    if test $status -eq 0; then
	FITSIO_LIB_FLAG=$lib_flag
	VERSION=$(./tst.exe)
	break
    fi
done
rm -f tst.c tst.exe

if test "$FITSIO_LIB_FLAG" == "?"; then
    cat <<EOF

I cannot find the fitsio library.  If you have it, modify the Makefile
to reflect your system configuration

EOF
    FITSIO_INCLUDE_FLAG="-I/path_to_the_folder_that_contains_fitsio_h"
    FITSIO_LIB_FLAG="-L/path_to_the_folder_that_contains_libcfitsio"
else
    echo ""
    echo "Found fitsio library version " $VERSION
    echo ""
fi


sed "s!%%FITSIO_INCLUDE_FLAG%%!${FITSIO_INCLUDE_FLAG}!g" Makefile.in > tmp
sed "s!%%FITSIO_LIB_FLAG%%!${FITSIO_LIB_FLAG}!g" tmp > Makefile

rm -f tmp

exit 0