#!/bin/bash

# Generate a nice table of all the functional groups 
# from the html output of align functional groups
# just run in the same directory as the make groups
# and move the svg to an img/subdirectory

COUNTER=0

rm functional_groups.html

OUTFILE=out.tmp

cat << EOF > $OUTFILE

<!DOCTYPE html>

<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta name="description" content="List of functional groups">
    <meta name="author" content="Thomas D Daff">

    <title>Functional group index</title>

    <!-- Bootstrap core CSS -->
    <link rel="stylesheet"
          href="http://netdna.bootstrapcdn.com/bootstrap/3.0.0/css/bootstrap.min.css"
          type="text/css">

</head>

<body>
<!-- Wrap all page content here -->

<div id="wrap">
    <!-- Begin page content -->

    <div class="container">
        <div class="page-header">
            <h1>Functional group list</h1>
        </div>

        <p class="lead">All currently available groups are included in the
            table</p>
    </div>

    <div class="container">

        <table id="grouptable"
               class="table table-bordered table-hover table-condensed">
            <thead>
            <tr>
                <th>Identifier</th>
                <th>Information</th>
                <th>Image</th>
                <th></th>
                <th>Identifier</th>
                <th>Information</th>
                <th>Image</th>
            </tr>
            </thead>
            <tbody>

            <!-- START GROUPS -->
EOF


for html in `ls -tr *.html`
do
    COUNTER=$((COUNTER + 1))
    if [[ $(($COUNTER % 2)) == 1 ]]
    then
        echo "            <tr>" >> $OUTFILE
        cat $html >> $OUTFILE
        echo "                <td></td>" >> $OUTFILE
    else
        cat $html >> $OUTFILE
        echo "            </tr>" >> $OUTFILE
    fi
done

# close table row if necessary

if [[ $(($COUNTER % 2)) == 1 ]]
then
    echo "            </tr>" >> $OUTFILE
fi


cat << EOF >> $OUTFILE
            <!-- END GROUPS -->
            </tbody>
        </table>


    </div>

</div>

</body>
</html>
EOF

mv $OUTFILE functional_groups.html

# Put the required images in the subdir
mkdir img
mv *.svg img/


