name := "ranola"

version := "0.1.0"

scalaVersion := "2.11.7"

licenses += ("Apache-2.0", url("http://opensource.org/licenses/Apache-2.0"))

libraryDependencies ++= Seq(
  "org.scalanlp" %% "breeze" % "0.11.2",
  "org.scalanlp" %% "breeze-natives" % "0.11.2",
  "org.scalatest" % "scalatest_2.11" % "2.2.4"
)

resolvers ++= Seq(
  "Sonatype Snapshots" at "https://oss.sonatype.org/content/repositories/snapshots/",
  "Sonatype Releases" at "https://oss.sonatype.org/content/repositories/releases/"
)
