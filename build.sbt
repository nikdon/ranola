name := "ranola"

version := "1.0"

scalaVersion := "2.11.5"

libraryDependencies ++= Seq(
  "org.scalanlp" % "breeze_2.11" % "0.11.1",
  "org.scalanlp" % "breeze-natives_2.11" % "0.11.1",
  "org.scalatest" % "scalatest_2.11" % "2.2.4"
)

resolvers ++= Seq(
  // other resolvers here
  // if you want to use snapshot builds (currently 0.11-SNAPSHOT), use this.
  "Sonatype Snapshots" at "https://oss.sonatype.org/content/repositories/snapshots/",
  "Sonatype Releases" at "https://oss.sonatype.org/content/repositories/releases/"
)
