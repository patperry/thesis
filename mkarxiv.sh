#!/bin/bash 

ARXIV=arXiv

TITLE=perry-thesis



CHAPTERS=( \
    introduction \
    multivariate \
    lowrank \
    intrinsic-rank \
    cv-unsupervised \
    bcv-theory \
    projections \
    weighted-sums \
)

AUXFILES=( \
    abstract.tex \
    acknowledgment.tex \
    bibliography.tex \
    perry-thesis.bbl \
    perry-thesis.tex \
    perry-thesis-macros.sty \
    suthesis-2e.sty
)

DISTDIR=${ARXIV}/${TITLE}


mkdir -p "${DISTDIR}"

for CHAPTER in "${CHAPTERS[@]}"
do
    DSTDIR="${DISTDIR}/${CHAPTER}"

    mkdir -p ${DSTDIR}
    cp ${CHAPTER}/${CHAPTER}.tex ${DSTDIR}

    if [ -e "${CHAPTER}/plots" ]; then
        cp -r ${CHAPTER}/plots ${DSTDIR}
    fi
done

for AUXFILE in "${AUXFILES[@]}"
do
    cp ${AUXFILE} ${DISTDIR}
done 

pushd ${ARXIV} > /dev/null
rm -f ${TITLE}.tar.gz
tar czvf ${TITLE}.tar.gz . > /dev/null
popd > /dev/null

