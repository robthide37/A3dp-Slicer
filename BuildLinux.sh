#!/bin/bash

export ROOT=`pwd`
export NCORES=`nproc`
FOUND_GTK2=$(dpkg -l libgtk* | grep gtk2)
FOUND_GTK3=$(dpkg -l libgtk* | grep gtk-3)

unset name
while getopts ":hwdrigbsyu" opt; do
  case ${opt} in
    u )
        UPDATE_LIB="1"
        ;;
    i )
        BUILD_IMAGE="1"
        ;;
    d )
        BUILD_DEPS="1"
        ;;
    s )
        BUILD_SLIC3R="1"
        ;;
    t)
        BUILD_TESTS="1"
        ;;
    b )
        BUILD_DEBUG="1"
        ;;
    g )
        FOUND_GTK3=""
        ;;
    w )
        BUILD_WIPE="1"
        ;;
    r )
        BUILD_CLEANDEPEND="1"
        ;;
    h ) echo "Usage: ./BuildLinux.sh [-h][-w][-d][-r][-i][-g][-b][-s][-t][-u]"
        echo "   -h: this message"
        echo "   -w: wipe build directories before building"
        echo "   -d: build deps (optional)"
        echo "   -r: clean dependencies building files (reduce disk usage)"
        echo "   -i: Generate appimage (optional)"
        echo "   -g: force gtk2 build"
        echo "   -b: build with debug symbols"
        echo "   -s: build Slic3r/SuperSlicer"
        echo "   -t: build tests (in combination with -s)"
        echo "   -u: only update clock & dependency packets (optional and need sudo)"
        echo "For a first use, you want to 'sudo ./BuildLinux.sh -u'"
        echo "   and then './BuildLinux.sh -dsi'"
        exit 0
        ;;
  esac
done

if [ $OPTIND -eq 1 ]
then
    echo "Usage: ./BuildLinux.sh [-h][-w][-d][-r][-i][-g][-b][-s][-t][-u]"
    echo "   -h: this message"
    echo "   -w: wipe build directories before building"
    echo "   -d: build deps (optional)"
    echo "   -r: clean dependencies building files (reduce disk usage)"
    echo "   -i: Generate appimage (optional)"
    echo "   -g: force gtk2 build"
    echo "   -b: build with debug symbols"
    echo "   -s: build Slic3r/SuperSlicer"
    echo "   -t: build tests (in combination with -s)"
    echo "   -u: only update clock & dependency packets (optional and need sudo)"
    echo "For a first use, you want to 'sudo ./BuildLinux.sh -u'"
    echo "   and then './BuildLinux.sh -dsi'"
    exit 0
fi

if [[ -n "$FOUND_GTK3" ]]
then
    echo "Found GTK3"
else
    if [[ -n "$FOUND_GTK2" ]]
    then
        echo "Found GTK2"
    fi
fi

if [[ -n "$UPDATE_LIB" ]]
then
    echo -n -e "Updating linux ...\n"
    hwclock -s
    apt update
	apt install g++ m4
    if [[ -z "$FOUND_GTK3" ]]
    then
        echo -e "\nInstalling: libgtk2.0-dev libglew-dev libudev-dev libdbus-1-dev cmake git gettext fuse\n"
        apt install libgtk2.0-dev libglew-dev libudev-dev libdbus-1-dev cmake git gettext fuse
    else
        echo -e "\nFind libgtk-3, installing: libgtk-3-dev libglew-dev libudev-dev libdbus-1-dev cmake git gettext fuse\n"
        apt install libgtk-3-dev libglew-dev libudev-dev libdbus-1-dev cmake git gettext fuse
    fi
    # for ubuntu 22.04:
    ubu_version="$(cat /etc/issue)" 
    if [[ $ubu_version == "Ubuntu 22.04"* ]]
    then
        apt install curl libssl-dev libcurl4-openssl-dev m4
    fi
    if [[ -n "$BUILD_DEBUG" ]]
    then
        echo -e "\nInstalling: libssl-dev libcurl4-openssl-dev\n"
        apt install libssl-dev libcurl4-openssl-dev
    fi
    echo -e "done\n"
    exit 0
fi

FOUND_GTK2_DEV=$(dpkg -l libgtk* | grep gtk2.0-dev)
FOUND_GTK3_DEV=$(dpkg -l libgtk* | grep gtk-3-dev)
echo "FOUND_GTK2=$FOUND_GTK2)"
if [[ -z "$FOUND_GTK2_DEV" ]]
then
if [[ -z "$FOUND_GTK3_DEV" ]]
then
    echo "Error, you must install the dependencies before."
    echo "Use option -u with sudo"
    exit 0
fi
fi

echo "[1/9] Updating submodules..."
{
    # update submodule profiles
    pushd resources/profiles
    git submodule update --init
    popd
}

echo "[2/9] Changing date in version..."
{
    # change date in version
    sed -i "s/+UNKNOWN/_$(date '+%F')/" version.inc
}
echo "done"

# mkdir build
if [ ! -d "build" ]
then
    mkdir build
fi

# mkdir in deps
if [ ! -d "deps/build" ]
then
    mkdir deps/build
fi

if [[ -n "$BUILD_DEPS" ]]
then
    if [[ -n $BUILD_WIPE ]]
    then
       echo -e "\n wiping deps/build directory ...\n"
       rm -fr deps/build
       echo -e " ... done\n"
    fi
    # mkdir in deps
    if [ ! -d "deps/build" ]
    then
        mkdir deps/build
    fi
    echo "[3/9] Configuring dependencies..."
    BUILD_ARGS=""
    if [[ -n "$FOUND_GTK3_DEV" ]]
    then
        BUILD_ARGS="-DDEP_WX_GTK3=ON"
    fi
    if [[ -n "$BUILD_DEBUG" ]]
    then
        # have to build deps with debug & release or the cmake won't find evrything it needs
        mkdir deps/build/release
        pushd deps/build/release
            cmake ../.. -DDESTDIR="../destdir" $BUILD_ARGS
            make -j$NCORES
        popd
        BUILD_ARGS="${BUILD_ARGS} -DCMAKE_BUILD_TYPE=Debug"
    fi
    
    # cmake deps
    pushd deps/build
        cmake .. $BUILD_ARGS
        echo "done"
        
        # make deps
        echo "[4/9] Building dependencies..."
        make -j$NCORES
        echo "done"
        
        # rename wxscintilla
        echo "[5/9] Renaming wxscintilla library..."
        pushd destdir/usr/local/lib
            if [[ -z "$FOUND_GTK3_DEV" ]]
            then
                cp libwxscintilla-3.2.a libwx_gtk2u_scintilla-3.2.a
            else
                cp libwxscintilla-3.2.a libwx_gtk3u_scintilla-3.2.a
            fi
        popd
        echo "done"
        
    popd
    echo "done"
fi

# clean deps
if [[ -n "$BUILD_CLEANDEPEND" ]]
then
    echo -e "[6/9] Cleaning dependencies...\n"
    pushd deps/build
    pwd
    rm -fr dep_*
    popd > /dev/null
    echo -e "\n ... done\n"
fi

if [[ -n "$BUILD_SLIC3R" ]]
then
    echo "[7/9] Configuring Slic3r..."
    BUILD_ARGS=""
    if [[ -n "$FOUND_GTK3_DEV" ]]
    then
        BUILD_ARGS="-DSLIC3R_GTK=3"
    fi
    if [[ -n "$BUILD_DEBUG" ]]
    then
        BUILD_ARGS="${BUILD_ARGS} -DCMAKE_BUILD_TYPE=Debug"
    fi
    if [[ -n "$BUILD_TESTS" ]]
    then
        BUILD_ARGS="${BUILD_ARGS} -DCMAKE_BUILD_TESTS=1"
    else
        BUILD_ARGS="${BUILD_ARGS} -DCMAKE_BUILD_TESTS=0"
    fi
    
    # cmake
    pushd build
        cmake .. -DCMAKE_PREFIX_PATH="$PWD/../deps/build/destdir/usr/local" -DSLIC3R_STATIC=1 ${BUILD_ARGS}
        echo "done"
        
        #make avrdude-slic3r
        make avrdude-slic3r
        
        # make Slic3r
        echo "[8/9] Building Slic3r..."
        make -j$NCORES Slic3r

        # make .mo
        make gettext_po_to_mo
        
        # make OCCTWrapper.so
        make OCCTWrapper
        
        # update the pot
        make gettext_make_pot
    popd
    echo "done"
fi

# Give proper permissions to script
chmod 755 $ROOT/build/src/BuildLinuxImage.sh

echo "[9/9] Generating Linux app..."
    pushd build
        if [[ -n "$BUILD_IMAGE" ]]
        then
            $ROOT/build/src/BuildLinuxImage.sh -i
        else
            $ROOT/build/src/BuildLinuxImage.sh
        fi
    popd
echo "done"
