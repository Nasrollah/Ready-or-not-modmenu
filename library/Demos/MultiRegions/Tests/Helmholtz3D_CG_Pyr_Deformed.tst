<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Helmholtz 3D CG for Pyramid (deformed)</description>
    <executable>Helmholtz3D</executable>
    <parameters>Helmholtz3D_Pyr_Deformed.xml</parameters>
    <files>
        <file description="Session File">Helmholtz3D_Pyr_Deformed.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">2.76145e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">0.000161114</value>
        </metric>
    </metrics>
</test>
