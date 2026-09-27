# Build the existing Jekyll site, including its remote theme, into excluded test output.
require 'jekyll'
require 'jekyll-remote-theme'
require 'jekyll-seo-tag'
require 'jekyll-include-cache'
Dir.chdir(File.expand_path('..', __dir__))
config = Jekyll.configuration(
  'source' => Dir.pwd,
  'destination' => File.join(Dir.pwd, 'tests/.pcr-output/site'),
  'incremental' => false,
  'sass' => { 'cache_location' => File.join(Dir.pwd, 'tests/.pcr-tools/sass-cache') }
)
site = Jekyll::Site.new(config)
site.process
puts "Jekyll #{Jekyll::VERSION} built #{site.pages.length} pages to #{site.dest}"
