if [[ $(hostname) = *eslogin* ]]; then
  echo "Hazelhen detected, calling hazelhen build"
  . hazelhen_build.sh
else
  ccc && cmake -DOPENCMISS_BUILD_TYPE=RELEASE -DCMAKE_BUILD_TYPE=RELEASE  .. && make clean && make all
fi

