<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>SVV Quad P=12</description>
    <executable>ADRSolver</executable>
    <parameters>SVV_Quad.xml</parameters>
    <files>
        <file description="Session File">SVV_Quad.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="0">
            <value variable="u" tolerance="1e-05">0.0337347</value>
        </metric>
        <metric type="Linf" id="1">
            <value variable="u" tolerance="1e-05"> 0.205011</value>
        </metric>
    </metrics>
</test>
