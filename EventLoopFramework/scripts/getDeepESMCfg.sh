#!/bin/bash

GITHUB_URL=https://github.com/StealthStop
REPO_NAME=DeepESMCfg
RELEASE_URL="$GITHUB_URL/$REPO_NAME/releases"

STARTING_DIR=$PWD
CFG_DIRECTORY=$PWD
CHECKOUT_DIRECTORY=$PWD
TAG=
NO_SOFTLINK=
OVERWRITE=
VERBOSE=

# A POSIX variable
OPTIND=1    # Reset in case getopts has been used previously in the shell.

SOFTLINK_NAME=DeepEventShape
SOFTLINK_NAME_2=DeepEventShape_NonIsoMuon
SOFTLINK_SUFFIX=

function print_help {
    echo ""
    echo "Usage:"
    echo "    getTaggerCfg.sh -t RELEASE_TAG [-d checkout_directory] [-f cfg_filename] [-o] [-n] [-v]"
    echo ""
    echo "Options:"
    echo "    -t RELEASE_TAG :         This is the github release tag to check out (required option)"
    echo "    -d checkout_directory :  This is the directory where the configuration files will be downloaded to (default: .)"
    echo "    -f cfg_filename :        Specify this option to name the softlink to the cfg file something other than \"TopTagger.cfg\""
    echo "    -o :                     Overwrite the softlinks if they already exist"
    echo "    -l checkout location :   Location to check out tagger cfg files (default: .)"
    echo "    -n :                     Download files without producing softlinks"
    echo "    -v :                     increase verbosity: print more stuff... for those who like stuff"
    echo "    -s softlink_suffix :     Can add a suffix to the softlink name"
    echo ""
    echo "Description:"
    echo "    This script automatically downloads the top tagger configuration file and MVA training files (if necessary)"
    echo "    and produces a softlink to this file in your corrent directory.  This script should be run from the directory where"
    echo "    the tagger code will be run from.  Tagger configuration releases can be browsed at"
    echo "    $RELEASE_URL"
    echo ""
}

function print_ok {
    echo "  ______    __  ___" 
    echo " /  __  \\  |  |/  /" 
    echo "|  |  |  | |  '  / " 
    echo "|  |  |  | |    <  " 
    echo "|  \`--'  | |  .  \\ " 
    echo " \\______/  |__|\\__\\" 
    echo ""
}

# Initialize our own variables:

while getopts "h?d:f:F:t:l:s:nov" opt; do
    case "$opt" in
    h|\?)
        print_help
        exit 0
        ;;
    d)  CFG_DIRECTORY=$OPTARG
        ;;
    f)  SOFTLINK_NAME=$OPTARG
        ;;
    F)  SOFTLINK_NAME_2=$OPTARG
        ;;
    t)  TAG=$OPTARG
        ;;
    l)  CHECKOUT_DIRECTORY=$OPTARG
        ;;
    s)  SOFTLINK_SUFFIX=$OPTARG
        ;;
    n)  NO_SOFTLINK=NO
        ;;
    o)  OVERWRITE="-f"
        ;;
    v)  VERBOSE=1
        ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

if [[ -z $TAG ]]
then
    print_help
    exit 0
fi

if [[ -z $SOFTLINK_SUFFIX ]]
then
    SOFTLINK_NAME="${SOFTLINK_NAME}.cfg"
    SOFTLINK_NAME_2="${SOFTLINK_NAME_2}.cfg"
else
    SOFTLINK_NAME="${SOFTLINK_NAME}_${SOFTLINK_SUFFIX}.cfg"
    SOFTLINK_NAME_2="${SOFTLINK_NAME_2}_${SOFTLINK_SUFFIX}.cfg"
fi

echo " - Running getDeepESMCfg.sh"

# get source directory of bash script
# used for "Easter Egg"...
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

# check if OVERWRITE is set
# if OVERWRITE is set, ask user for confirmation before continuing
if [[ -z $OVERWRITE ]]
then
    # OVERWRITE is not set
    # check if softlink exists
    if [[ -L $SOFTLINK_NAME ]]
    then
        echo "INFO: OVERWRITE is not set. Existing softlinks will not be replaced."
        printf "  Enter \"o\" to overwrite existing softlinks, and anything else to continue without overwriting: "
        read answer
        if [[ $answer == "ok" ]]
        then
            # "Easter Egg"...
            print_ok
        elif [[ $answer == "o" ]]
        then
            OVERWRITE="-f"
        fi
    fi
else
    # OVERWRITE is set
    # check if file exists and is not a softlink
    if [[ -f $SOFTLINK_NAME && ! -L $SOFTLINK_NAME ]]
    then
        # ask user for confirmation before continuing
        echo   "INFO: OVERWRITE is set. Existing files will be replaced."
        printf "  Enter (Y/y/yes/si/oui/ja/da) to replace existing files, and anything else to quit: "
        read answer
        if [[ $answer == "ok" ]]
        then
            # "Easter Egg"...
            print_ok
            exit 0
        fi
        if [[ $answer == "Y" || $answer == "y" || $answer == "yes" || $answer == "si" || $answer == "oui" || $answer == "ja" || $answer == "da" ]]
        then
            echo " - Continuing..."
        else
            echo " - Quitting..."
            exit 0
        fi
    fi
fi


# Check that CFG_DIRECTORY is a directory
if [ ! -d $CFG_DIRECTORY ]
then
    echo $CFG_DIRECTORY " Is not a valid directory!"
    exit 1
fi

cd $CFG_DIRECTORY

if [ ! -d ${REPO_NAME}_$TAG ]
then
    echo " - Downloading this REPO-TAG: ${REPO_NAME}_$TAG"
    if [[ ! -z $VERBOSE ]] # True if VERBOSE is set
    then
        wget $GITHUB_URL/$REPO_NAME/archive/$TAG.tar.gz
    else
        wget $GITHUB_URL/$REPO_NAME/archive/$TAG.tar.gz &> /dev/null
    fi
    if [ -f $TAG.tar.gz ]
    then
        tar xzf $TAG.tar.gz
        rm $TAG.tar.gz
        mv $REPO_NAME-$TAG ${REPO_NAME}_$TAG
    else
        echo "ERROR: Failed to download $GITHUB_URL/$REPO_NAME/archive/$TAG.tar.gz"
        echo "  Check that the REPO-TAG that you entered (${REPO_NAME}_$TAG)"
        echo "  exists at $RELEASE_URL"
        echo "  Check your spelling... you may have a typo! Copy and paste are your friends."
        exit 1
    fi
else
    echo " - Skipping the download of the requested REPO-TAG because the directory "${REPO_NAME}_$TAG" is already present"
fi


cd ${REPO_NAME}_$TAG
DOWNLOAD_DIR=$PWD

if [[ ! -z $VERBOSE ]] # True if VERBOSE is set
then
    echo "INFO: DOWNLOAD_DIR is $DOWNLOAD_DIR"
fi

MVAFILES=

if [ -f DeepEventShape.cfg ]
then
    MVAFILES=$(grep "modelFile" DeepEventShape.cfg | sed 's/[^"]*"\([^"]*\)"/\1/')
    MISSING=
    if [[ ! -z ${MVAFILES// } ]]
    then
        for MVAFILE in $MVAFILES; do
            if [ ! -f $MVAFILE ]
            then
                MISSING="yes"
                break
            fi
        done
        if [[ ! -z ${MISSING// } ]]
        then
            echo " - Downloading MVA files"
            MVATARBALL=MVAFILES.tar.gz
            if [[ ! -z $VERBOSE ]] # True if VERBOSE is set
            then
                wget $GITHUB_URL/$REPO_NAME/releases/download/$TAG/$MVATARBALL
            else
                wget $GITHUB_URL/$REPO_NAME/releases/download/$TAG/$MVATARBALL &> /dev/null
            fi
            if [ ! -f $MVATARBALL ]
            then
                echo "ERROR: MVA tarball "$MVATARBALL" not found!!!"
                MVATARBALL=${MVAFILES%.*}.tar.gz
                echo "Now trying $MVATARBALL instead"
                if [[ ! -z $VERBOSE ]] # True if VERBOSE is set
                then
                    wget $GITHUB_URL/$REPO_NAME/releases/download/$TAG/$MVATARBALL
                else
                    wget $GITHUB_URL/$REPO_NAME/releases/download/$TAG/$MVATARBALL &> /dev/null
                fi
                if [ ! -f $MVATARBALL ]
                then
                    echo "ERROR: MVA tarball "$MVATARBALL" not found!!!"
                    exit 1
                fi
            fi
            tar xzf $MVATARBALL
            rm $MVATARBALL
        fi
    fi
fi

# make all files in DOWNLOAD_DIR read only
# a (all) = ugo (user group others)
chmod a-w *

cd $STARTING_DIR

if [[ ! -z $VERBOSE ]] # True if VERBOSE is set
then
    echo "INFO: STARTING_DIR is $STARTING_DIR"
fi

# If OVERWRITE is set, make solftlinks (using ln) with -f
# If OVERWRITE is not set, make solftlinks (using ln)
# Pipe output to /dev/null

# Note: "> /dev/null 2>&1" does this:
# stdin  ==> fd 0      (default fd 0)
# stdout ==> /dev/null (default fd 1)
# stderr ==> stdout    (default fd 2)

# [[ -z STRING ]] : True if the length of "STRING" is zero, False if "STRING" has nonzero length
if [[ -z $NO_SOFTLINK ]]
then
    # create softlinks
    ln $OVERWRITE -s $DOWNLOAD_DIR/DeepEventShape.cfg $CHECKOUT_DIRECTORY/$SOFTLINK_NAME > /dev/null 2>&1 && echo " - Created softlinks to $REPO_NAME config file"
    ln $OVERWRITE -s $DOWNLOAD_DIR/DeepEventShape_NonIsoMuon.cfg $CHECKOUT_DIRECTORY/$SOFTLINK_NAME_2 > /dev/null 2>&1 && echo " - Created softlinks to $REPO_NAME config file"
    if [[ ! -z ${MVAFILES// } ]] 
    then
        for MVAFILE in $MVAFILES; do
            MVAFILE_SOFTLINK_NAME=$MVAFILE
            if [[ ! -z $SOFTLINK_SUFFIX ]]
            then
                IFS='.' read -ra nameList <<< "$MVAFILE"
                MVAFILE_SOFTLINK_NAME="${nameList[0]}_${SOFTLINK_SUFFIX}.${nameList[1]}"
            fi
            ln $OVERWRITE -s $DOWNLOAD_DIR/$MVAFILE $CHECKOUT_DIRECTORY/$MVAFILE_SOFTLINK_NAME > /dev/null 2>&1 && echo " - Created softlinks to $REPO_NAME MVA files"
        done
    fi
fi
