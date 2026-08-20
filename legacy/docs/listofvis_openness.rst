.. _listofvis_openness:

Openness
========

Positive openness
-----------------

.. image:: ./figures/rvtvis_qgis_Pivola_dem_05m_OPEN-POS_R10_D16_8bit.jpg
   :width: 650px

Openness is another proxy for relief shading (see also: sky-view factor). The method is based on estimating the mean horizon elevation angle within a defined search radius. The mean value of all zenith angles gives positive openness, while the mean nadir value gives negative openness. 

Openness considers the whole sphere for calculation, not just the celestial hemisphere as SVF does. The result of this is a much ‘flatter’ image, devoid of general topography—a kind of trend-removed image. As the visual impression of the general topography is lost, interpretation becomes a bit trickier. 

However, openness has big advantage for automatic feature detection because ‘signatures’ of features are more homogeneous because they are the same irrespectively of their location on a plane or slope.

Negative openness
-----------------

.. image:: ./figures/rvtvis_qgis_Pivola_dem_05m_OPEN-NEG_R10_D16_8bit.jpg
   :width: 650px

Negative openness is not the inverse of positive openness and it provides additional information. While positive openness highlights topographic convexities (e.g. ridges between hollow ways and rims of bomb craters), negative openness emphasizes the lowest parts of concavities, (e.g. the actual hollow ways, the lowest parts of gorges and the lower edges of cliffs). 

For consistent readability, it is recommended that negative openness is displayed with inverted greyscale (i.e. darker for higher values), so that concave features are always presented by dark tones.

Positive and negative openness are very useful to highlight positive and negative relief features, respectively. As openness removes the visual impression of overall landscape forms, it is not affected by saturation due to gentle or steep slopes and may be used in a varied topography. Because of the ability to differentially highlight positive and negative relief features, it is particularly suitable for targeted detection of these features.

.. seealso:: :ref:`listofvis_svf`.
