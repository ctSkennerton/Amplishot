#Amplishot



##Install

From github run 
```
git clone git://github.com/ctSkennerton/Amplishot.git
cd Amplishot
git submodule init
git submodule update
python setup.py install
```

You will need to add the installed location 
i.e. $PREFIX/lib/python2.7/site-packages/ into your LD_LIBRARY_PATH as
Ampliphot includes c extensions

```
export LD_LIBRARY_PATH="$PREFIX/lib/python2.7/site-packages/:$LD_LIBRARY_PATH"
```

to make the changes permenent add the above line into your bashrc
