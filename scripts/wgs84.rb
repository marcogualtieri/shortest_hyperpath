#!/usr/bin/env ruby
#
# Various procedures to handle geographical coordinates.
#
# Code taken from various sources in other languages and converted
# into ruby.
#
# Copyright:: Copyright (c) 2012 Mirko Maischberger
# License:: Distributed under the same terms as Ruby
  
module Wgs84
  A0 = 6367449.14580084
  B0 = 16038.4295531591
  C0 = 16.8326133343344
  D0 = 0.0219844042737573
  E0 = 0.0003127051795045
  ELL_A = 6378137.0
  ELL_B = 6356752.314
  UTM_SCALE = 0.9996
  ECC_SQUARED = 0.00669437999013
end

# convert degrees (decimal format) to radiants 
def rad(deg)
  return (deg.to_f*Math::PI)/180.0
end

# Calculate geodesic distance (in m) between two points specified by
# latitude/longitude (in numeric degrees) using Vincenty inverse
# formula for ellipsoids.
#
# Returns: Vicenti distance from lat1, lon1 to lat2, lon2 on the
# WGS-84 ellipsoid.
#
def vicenti_distance(lat1, lon1, lat2, lon2)
  a,b,f = 6378137.0, 6356752.3142, 1.0/298.257223563 # WGS-84 ellipsiod
  l = rad(lon2-lon1)
  u1 = Math::atan((1.0-f) * Math::tan(rad(lat1)))
  u2 = Math::atan((1.0-f) * Math::tan(rad(lat2)))
  sinU1, cosU1 = Math::sin(u1), Math::cos(u1)
  sinU2, cosU2 = Math::sin(u2), Math::cos(u2)
  
  lambda, lambdaP = l, 2.0*Math::PI
  iterLimit = 20;
  while ((lambda-lambdaP).abs > 1e-12 && --iterLimit>0) do 
    sinLambda, cosLambda = Math::sin(lambda), Math::cos(lambda)
    sinSigma = Math::sqrt((cosU2*sinLambda) * (cosU2*sinLambda) + 
                          (cosU1*sinU2-sinU1*cosU2*cosLambda) * 
                          (cosU1*sinU2-sinU1*cosU2*cosLambda))
    if (sinSigma==0) 
      return 0.0  # co-incident points
    end
    cosSigma = sinU1*sinU2 + cosU1*cosU2*cosLambda
    sigma = Math::atan2(sinSigma, cosSigma)
    sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma
    cosSqAlpha = 1.0 - sinAlpha*sinAlpha
    cos2SigmaM = cosSqAlpha == 0.0 ? 0.0 : cosSigma - 2.0*sinU1*sinU2/cosSqAlpha
    c = f/16.0*cosSqAlpha*(4.0+f*(4.0-3.0*cosSqAlpha))
    lambdaP = lambda
    lambda = (l + (1.0-c) * f * sinAlpha *
              (sigma + c*sinSigma*(cos2SigmaM+c*cosSigma*(-1.0+2.0*cos2SigmaM*cos2SigmaM))))
  end
  
  # formula failed to converge
  if (iterLimit == 0) 
    return nil
  end
  uSq = cosSqAlpha * (a*a - b*b) / (b*b)
  ag = 1.0 + uSq/16384.0*(4096.0+uSq*(-768.0+uSq*(320.0-175.0*uSq)))
  bg = uSq/1024.0 * (256.0+uSq*(-128.0+uSq*(74.0-47.0*uSq)))
  deltaSigma = bg*sinSigma*(cos2SigmaM+bg/4.0*
                            (cosSigma*(-1.0+2.0*cos2SigmaM*cos2SigmaM)-
                             bg/6.0*cos2SigmaM*(-3.0+4.0*sinSigma*sinSigma)*
                             (-3.0+4.0*cos2SigmaM*cos2SigmaM)))
  return b*ag*(sigma-deltaSigma)
end

## WGS2UTM taken from JS version by from Charles L. Taylor
## http://home.hiwaay.net/~taylorc/toolbox/geography/geoutm.html

class Utm
  attr_accessor :northing, :easting, :zone
end

# Computes the ellipsoidal distance from the equator to a point at a
# given latitude.
#
# Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
# GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994.
#
# phi - Latitude of the point, in radians.
#
def meridian_arc_len(phi)
  n = (Wgs84::ELL_A - Wgs84::ELL_B) / (Wgs84::ELL_A + Wgs84::ELL_B);
  alpha = ((Wgs84::ELL_A + Wgs84::ELL_B) / 2.0) *
    (1.0 + ((n**2.0) / 4.0) + ((n**4.0) / 64.0));
  
  beta = (-3.0 * n / 2.0) + (9.0 * (n**3.0) / 16.0) +
    (-3.0 * (n**5.0) / 32.0);
  gamma = (15.0 * (n**2.0) / 16.0) +
    (-15.0 * (n**4.0) / 32.0);
  delta = (-35.0 * (n**3.0) / 48.0) +
    (105.0 * (n**5.0) / 256.0);
  
  epsilon = (315.0 * (n**4.0) / 512.0);
  
  alpha * (phi + (beta * Math.sin(2.0 * phi)) +
           (gamma * Math.sin(4.0 * phi)) +
           (delta * Math.sin(6.0 * phi)) +
           (epsilon * Math.sin(8.0 * phi)));
end

# Central meridian of specified UTM zone
def central_meridian(zone)
  rad(-183.0 + (zone * 6.0))
end

# Converts a radiants latitude/longitude pair to x and y coordinates in the
# Transverse Mercator projection.  Note that Transverse Mercator is not
# the same as UTM; a scale factor is required to convert between them.
# 
# Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
# GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994.
#
def wgs2xy(phi, lambda, zone)
  lambda0 = central_meridian(zone)

  ep2 = ((Wgs84::ELL_A**2.0) - (Wgs84::ELL_B**2.0)) /
    (Wgs84::ELL_B**2.0)
    
  nu2 = ep2 * (Math.cos(phi) ** 2.0)
  
  bign = (Wgs84::ELL_A**2.0) / (Wgs84::ELL_B*Math.sqrt(1 + nu2))
  
  t = Math.tan(phi)
  t2 = t*t
  tmp = (t2*t2*t2)-(t**6.0)
  
  l = lambda - lambda0
    
  l3coef = 1.0 - t2 + nu2
  l4coef = 5.0 - t2 + 9 * nu2 + 4.0 * (nu2**2.0)
  l5coef = 5.0 - 18.0*t2 + (t2**2.0) + 14.0*nu2 - 58.0*t2*nu2
  l6coef = 61.0 - 58.0*t2 + (t2**2.0) + 270.0*nu2 - 330.0*t2*nu2
  l7coef = 61.0 - 479.0 * t2 + 179.0 * (t2**2.0) - (t2**3.0)
  l8coef = 1385.0 - 3111.0*t2 + 543.0 * (t2**2.0) - (t2**3.0)

  utm = Utm.new

  utm.easting = bign * Math.cos(phi) * l +
    (bign / 6.0 * (Math.cos(phi)**3.0) * l3coef * (l**3.0)) +
    (bign / 120.0 * (Math.cos(phi)**5.0) * l5coef * (l**5.0)) +
    (bign / 5040.0 * (Math.cos(phi)**7.0) * l7coef * (l**7.0))
    
  utm.northing = meridian_arc_len(phi) +
    (t / 2.0 * bign * (Math.cos(phi)**2.0) * (l**2.0)) +
    (t / 24.0 * bign * (Math.cos(phi)**4.0) * l4coef * (l**4.0)) +
    (t / 720.0 * bign * (Math.cos(phi)**6.0) * l6coef * (l**6.0)) +
    (t / 40320.0 * bign * (Math.cos(phi)**8.0) * l8coef * (l**8.0))

  utm.zone = zone

  utm
end

#  Converts a latitude/longitude pair to x and y coordinates in the
#  Universal Transverse Mercator projection.
#
#  lat - Latitude of the point, in decimal degrees.
#  lon - Longitude of the point, in decimal degrees.
#
def wgs2utm(lat, lon)
  utm = wgs2xy(rad(lat), rad(lon), ((lon + 180.0) / 6.0).floor + 1)
  utm.easting = utm.easting * Wgs84::UTM_SCALE + 500000.0;
  utm.northing = utm.northing * Wgs84::UTM_SCALE;
  if utm.northing  < 0.0
    utm.northing  = utm.northing + 10000000.0;
  end
  utm
end

if __FILE__ == $0 then
  
  if ARGV.size == 4
    lat1, lon1, lat2, lon2 = ARGV.map { |a| a.to_f }
    puts "#{vicenti_distance(lat1, lon1, lat2, lon2)}m"
  elsif ARGV.size == 2
    lat, lon = ARGV.map { |a| a.to_f }
    utm = wgs2utm(lat, lon)
    puts "#{utm.northing};#{utm.easting};#{utm.zone}"
  else
    puts "use the source, luke!"
  end
  
end
